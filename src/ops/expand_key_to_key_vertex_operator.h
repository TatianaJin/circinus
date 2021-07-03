// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/expand_vertex_traverse_context.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

template <typename G, bool intersect_candidates = true>
class ExpandKeyToKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandKeyToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               SubgraphFilter* filter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, filter) {}

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, (ExpandVertexTraverseContext*)ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<ctx->type>(batch_size, (ExpandVertexTraverseContext*)ctx);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToKeyVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(uint32_t batch_size, ExpandVertexTraverseContext* ctx) const {
    auto data_graph = (G*)(ctx->current_data_graph);
    uint32_t output_num = 0;
    while (ctx->hasNextInput()) {
      std::vector<VertexID> new_keys;
      const auto& input = ctx->getCurrentInput();
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_[parents_[i]];
        DCHECK_LT(key, input.getNumKeys());
        uint32_t key_vid = input.getKeyVal(key);
        auto neighbors = data_graph->getOutNeighborsWithHint(key_vid, ALL_LABEL, i);
        if (i == 0) {
          if (!intersect_candidates) {
            removeExceptions(neighbors, &new_keys, exceptions);
          } else {
            intersect(*candidates_, neighbors, &new_keys, exceptions);
            if
              constexpr(isProfileMode(profile)) {
                updateIntersectInfo(candidates_->size() + neighbors.size(), new_keys.size());
              }
          }
        } else {
          auto new_keys_size = new_keys.size();
          (void)new_keys_size;
          intersectInplace(new_keys, neighbors, &new_keys);
          if
            constexpr(isProfileMode(profile)) {
              updateIntersectInfo(new_keys_size + neighbors.size(), new_keys.size());
            }
        }
        if (new_keys.size() == 0) {
          break;
        }
      }
      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
          if
            constexpr(isProfileWithMiniIntersectionMode(profile)) {
              // consider reuse of partial intersection results at each parent
              std::vector<VertexID> parent_tuple(parents_.size());
              for (uint32_t i = 0; i < parents_.size(); ++i) {
                uint32_t key_vid = input.getKeyVal(query_vertex_indices_[parents_[i]]);
                parent_tuple[i] = key_vid;
                ctx->updateDistinctSICount(i, parent_tuple, i);
              }
            }
        }
      if (new_keys.size() != 0) {
        for (VertexID new_key : new_keys) {
#ifdef USE_FILTER
          CompressedSubgraphs output(input, new_key, same_label_set_indices_, set_pruning_threshold_, false);
          if (output.empty() || filter(output)) continue;
#else
          CompressedSubgraphs output(input, new_key, same_label_set_indices_, set_pruning_threshold_);
          if (output.empty()) continue;
#endif
          ctx->outputs->emplace_back(std::move(output));
          ++output_num;
        }
      }
      ctx->nextInput();
      if (output_num >= batch_size) {
        break;
      }
      new_keys.clear();
    }
    return output_num;
  }
};

}  // namespace circinus

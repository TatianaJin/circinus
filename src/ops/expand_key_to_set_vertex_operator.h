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

#include <algorithm>
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

template <typename G>
class ExpandKeyToSetVertexOperator : public ExpandVertexOperator {
 public:
  ExpandKeyToSetVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               SubgraphFilter* filter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, filter) {}

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    if (ctx->query_type == QueryType::Profile) return expandInner<QueryType::Profile>(batch_size, ctx);
    CHECK(ctx->query_type == QueryType::ProfileWithMiniIntersection) << "Unknown query type "
                                                                     << (uint32_t)ctx->query_type;
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, ctx);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToSetVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first, input_key_size.second + 1};
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(uint32_t batch_size, TraverseContext* base_ctx) const {
    auto ctx = (ExpandVertexTraverseContext*)base_ctx;
    uint32_t output_num = 0;
    for (; output_num < batch_size && ctx->hasNextInput(); ctx->nextInput()) {
      const auto& input = ctx->getCurrentInput();
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      std::vector<VertexID> new_set;
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_.at(parents_[i]);
        uint32_t key_vid = input.getKeyVal(key);
        auto neighbors = ((G*)ctx->current_data_graph)->getOutNeighborsWithHint(key_vid, target_label_, i);
        if (i == 0) {
          intersect(*candidates_, neighbors, &new_set, exceptions);
          if
            constexpr(isProfileMode(profile)) {
              ctx->updateIntersectInfo(candidates_->size() + neighbors.size(), new_set.size());
            }
        } else {
          auto new_set_size = new_set.size();
          (void)new_set_size;
          intersectInplace(new_set, neighbors, &new_set);
          if
            constexpr(isProfileMode(profile)) {
              ctx->updateIntersectInfo(new_set_size + neighbors.size(), new_set.size());
            }
        }
        if (new_set.size() == 0) {
          break;
        }
      }
      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
          // consider reuse of partial intersection results at each parent
          if (isProfileWithMiniIntersectionMode(profile)) {
            std::vector<VertexID> parent_tuple(parents_.size());
            for (uint32_t j = 0; j < parents_.size(); ++j) {
              uint32_t key_vid = input.getKeyVal(query_vertex_indices_.at(parents_[j]));
              parent_tuple[j] = key_vid;
              ctx->updateDistinctSICount(j, parent_tuple, j);
            }
          }
        }
      if (!new_set.empty()) {
#ifdef USE_FILTER
        auto output = ctx->newOutput(input, std::move(new_set));
        if (filter(*output)) {
          ctx->popOutput();
          continue;
        }
#else
        auto output = ctx->newOutput(input, std::move(new_set), same_label_set_indices_, set_pruning_threshold_);
        if (output == nullptr) {
          continue;
        }
#endif
        ++output_num;
        // TODO(by) break if batch_size is reached
      }
    }
    return output_num;
  }
};

}  // namespace circinus

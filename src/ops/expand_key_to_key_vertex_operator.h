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
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandKeyToKeyVertexOperator : public ExpandVertexOperator {
  std::vector<unordered_set<std::string>> parent_tuple_sets_;

 public:
  ExpandKeyToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               SubgraphFilter* filter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, filter) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<QueryType::Execute>(outputs, batch_size);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<QueryType::Profile>(outputs, batch_size);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToKeyVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new ExpandKeyToKeyVertexOperator(*this);
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    if
      constexpr(isProfileMode(profile)) { parent_tuple_sets_.resize(parents_.size()); }
    uint32_t output_num = 0;
    while (input_index_ < current_inputs_->size()) {
      std::vector<VertexID> new_keys;
      const auto& input = (*current_inputs_)[input_index_];
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_[parents_[i]];
        DCHECK_LT(key, input.getNumKeys());
        uint32_t key_vid = input.getKeyVal(key);
        if(use_bipartite_graph_flag)current_data_graph_=bg_pointers_[i];// must only use validate things in BipartiteGraph then
        if (i == 0) {
          if(use_bipartite_graph_flag)removeExceptions(current_data_graph_->getOutNeighbors(key_vid), &new_keys, exceptions);
          else intersect(*candidates_, current_data_graph_->getOutNeighbors(key_vid), &new_keys, exceptions);
          if
            constexpr(isProfileMode(profile)) {
              updateIntersectInfo(candidates_->size() + current_data_graph_->getVertexOutDegree(key_vid),
                                  new_keys.size());
            }
        } else {
          auto new_keys_size = new_keys.size();
          intersectInplace(new_keys, current_data_graph_->getOutNeighbors(key_vid), &new_keys);
          if
            constexpr(isProfileMode(profile)) {
              updateIntersectInfo(new_keys_size + current_data_graph_->getVertexOutDegree(key_vid), new_keys.size());
            }
        }
        if (new_keys.size() == 0) {
          break;
        }
      }
      if
        constexpr(isProfileMode(profile)) {
          total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
          // consider reuse of partial intersection results at each parent
          std::vector<VertexID> parent_tuple(parents_.size());
          for (uint32_t i = 0; i < parents_.size(); ++i) {
            uint32_t key_vid = input.getKeyVal(query_vertex_indices_[parents_[i]]);
            parent_tuple[i] = key_vid;
            distinct_intersection_count_ +=
                parent_tuple_sets_[i].emplace((char*)parent_tuple.data(), (i + 1) * sizeof(VertexID)).second;
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
          outputs->emplace_back(std::move(output));
          ++output_num;
        }
      }
      ++input_index_;
      if (output_num >= batch_size) {
        break;
      }
      new_keys.clear();
    }
    return output_num;
  }
};

}  // namespace circinus

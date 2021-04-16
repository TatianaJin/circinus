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
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandKeyToSetVertexOperator : public ExpandVertexOperator {
  std::vector<unordered_set<std::string>> parent_tuple_sets_;

 public:
  ExpandKeyToSetVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               SubgraphFilter* filter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, filter) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<QueryType::Execute>(outputs, batch_size);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size,
                                 uint32_t query_type) override {
    if (query_type == 1) {
      return expandInner<QueryType::Profile>(outputs, batch_size);
    } else if (query_type == 2) {
      return expandInner<QueryType::ProfileWithMiniIntersection>(outputs, batch_size);
    }
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToSetVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new ExpandKeyToSetVertexOperator(*this);
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    for (; output_num < batch_size && input_index_ < current_inputs_->size(); ++input_index_) {
      const auto& input = (*current_inputs_)[input_index_];
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      std::vector<VertexID> new_set;
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_[parents_[i]];
        uint32_t key_vid = input.getKeyVal(key);
        if (i == 0) {
          intersect(*candidates_, current_data_graph_->getOutNeighbors(key_vid), &new_set, exceptions);
          if
            constexpr(isProfileMode(profile) || isProfileWithMiniIntersectionMode(profile)) {
              updateIntersectInfo(candidates_->size() + current_data_graph_->getVertexOutDegree(key_vid),
                                  new_set.size());
            }
        } else {
          auto new_set_size = new_set.size();
          intersectInplace(new_set, current_data_graph_->getOutNeighbors(key_vid), &new_set);
          if
            constexpr(isProfileMode(profile) || isProfileWithMiniIntersectionMode(profile)) {
              updateIntersectInfo(new_set_size + current_data_graph_->getVertexOutDegree(key_vid), new_set.size());
            }
        }
        if (new_set.size() == 0) {
          break;
        }
      }
      if
        constexpr(isProfileMode(profile) || isProfileWithMiniIntersectionMode(profile)) {
          total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
          // consider reuse of partial intersection results at each parent
          if (isProfileWithMiniIntersectionMode(profile)) {
            std::vector<VertexID> parent_tuple(parents_.size());
            parent_tuple_sets_.resize(parents_.size());
            for (uint32_t i = 0; i < parents_.size(); ++i) {
              uint32_t key_vid = input.getKeyVal(query_vertex_indices_[parents_[i]]);
              parent_tuple[i] = key_vid;
              distinct_intersection_count_ +=
                  parent_tuple_sets_[i].emplace((char*)parent_tuple.data(), (i + 1) * sizeof(VertexID)).second;
            }
          }
        }
      if (!new_set.empty()) {
#ifdef USE_FILTER
        CompressedSubgraphs output(input, std::move(new_set));
        if (filter(output)) {
          continue;
        }
#else
        CompressedSubgraphs output(input, std::move(new_set), same_label_set_indices_, set_pruning_threshold_);
        if (output.empty()) {
          continue;
        }
#endif
        outputs->emplace_back(std::move(output));
        ++output_num;
        // TODO(by) break if batch_size is reached
      }
    }
    return output_num;
  }
};

}  // namespace circinus

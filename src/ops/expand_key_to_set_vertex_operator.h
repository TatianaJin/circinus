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
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"

namespace circinus {

class ExpandKeyToSetVertexOperator : public ExpandVertexOperator {
 public:
  ExpandKeyToSetVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    while (input_index_ < current_inputs_->size()) {
      const auto& input = (*current_inputs_)[input_index_++];
      std::vector<VertexID> new_set;
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_[parents_[i]];
        uint32_t key_vid = input.getKeyVal(key);
        if (i == 0) {
          intersect(*candidates_, current_data_graph_->getOutNeighbors(key_vid), &new_set);
        } else {
          intersectInplace(new_set, current_data_graph_->getOutNeighbors(key_vid), &new_set);
        }
        if (new_set.size() == 0) {
          break;
        }
      }
      new_set.erase(std::remove_if(new_set.begin(), new_set.end(),
                                   [&input](VertexID set_vertex) { return input.isExisting(set_vertex); }),
                    new_set.end());
      if (new_set.size() != 0) {
        outputs->emplace_back(input, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
        ++output_num;
        // TODO(by) break if batch_size is reached
      }
    }
    return output_num;
  }
};

}  // namespace circinus

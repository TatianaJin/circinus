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

#include <vector>

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"

namespace circinus {

class ExpandKeyKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandKeyKeyVertexOperator(const QueryGraph* g, std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                       const std::vector<int>& cover_table, std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices) : ExpandVertexOperator(g, parents, target_vertex, cover_table, query_vertex_indices) {
  }

  void input(const std::vector<CompressedSubgraphs>* inputs, const Graph* data_graph) {
    inputs_ = inputs;
    data_graph_ = data_graph;
    inputs_idx_ = 0;
  }  

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    while (inputs_idx_ < inputs_->size()) {
      std::vector<VertexID> new_keys;
      const auto& input = (*inputs_)[inputs_idx_];
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_[parents_[i]];
        uint32_t key_vid = input.getKeyVal(key);
        if (i == 0) {
          intersect(*candidates_, data_graph_->getOutNeighbors(key_vid), &new_keys);
        } else {
          intersectInplace(new_keys, data_graph_->getOutNeighbors(key_vid), &new_keys);
        }
        if(new_keys.size() == 0) {
          break;
        }
      }
      if (new_keys.size() != 0) {
        for (VertexID new_key : new_keys) {
          outputs->emplace_back(input, new_key);
          output_num++;
        }
      }
      if (output_num >= batch_size) {
        break;
      }
      inputs_idx_++;
    }
    return output_num;
  }

 protected:
  const std::vector<CompressedSubgraphs>* inputs_;
  const Graph* data_graph_;
  uint32_t inputs_idx_;

};

}  // namespace circinus

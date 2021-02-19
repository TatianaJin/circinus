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
#include "ops/traverse_operator.h"

namespace circinus {

class ExpandIntoOperator : public TraverseOperator {
 public:
  ExpandIntoOperator(const QueryGraph* g, std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                       const std::vector<int>& cover_table, std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices) : g_(g), parents_(parents), target_vertex_(target_vertex), cover_table_(cover_table), query_vertex_indices_(query_vertex_indices) {
    // FIXME
  }

  void input(const std::vector<CompressedSubgraphs>* inputs, const Graph* data_graph) {
    inputs_ = inputs;
    data_graph_ = data_graph;
    inputs_idx_ = 0;
  }  
  
  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    while (inputs_idx_ < inputs_->size()) {
      std::unordered_map<VertexID, bool> visited;
      auto input = (*inputs_)[inputs_idx_];
      auto key_vertex_id = input.getKeyVal(query_vertex_indices_[target_vertex_]);
      auto key_out_neighbors = data_graph_->getOutNeighbors(key_vertex_id);
      CompressedSubgraphs new_output = input; 
      bool add = true;
      for (QueryVertexID vid : parents_) {
        std::vector<VertexID> new_set;
        uint32_t id = query_vertex_indices_[vid];
        intersect(*(input.getSet(id)), key_out_neighbors, &new_set);
        if (new_set.size() == 0) {
          add = false;
          break;
        }
        new_output.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
      }
      
      if (add) {
        outputs->emplace_back(new_output);
        output_num++;
      }
      inputs_idx_++;
      
      if (output_num >= batch_size) {
        break;
      }
    }
    return output_num;
  }
 
 protected:
  const QueryGraph* g_;
  const std::vector<int> cover_table_;
  std::unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;
  std::vector<QueryVertexID> parents_;
  QueryVertexID target_vertex_;
  const std::vector<CompressedSubgraphs>* inputs_;
  const Graph* data_graph_;
  uint32_t inputs_idx_;
};

}  // namespace circinus

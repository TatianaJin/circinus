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

class ExpandSetToKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandSetToKeyVertexOperator(const QueryGraph* g, std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
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
      std::unordered_map<VertexID, bool> visited;
      auto input = (*inputs_)[inputs_idx_];
      auto parent_match = input.getSet(query_vertex_indices_[parents_[0]]);
      for (VertexID vid : *parent_match.get()) {
        const auto& out_neighbors = data_graph_->getOutNeighbors(vid);
        for (uint32_t i = 0; i < out_neighbors.second; ++i) {
          VertexID key_vertex_id = *(out_neighbors.first + i);
          if (visited[key_vertex_id] == false || !isInCandidates(key_vertex_id)) {
            // TODO make sure key_vertex_id is in candidates_
            visited[key_vertex_id] = true;
            auto key_out_neighbors = data_graph_->getOutNeighbors(key_vertex_id);
            CompressedSubgraphs new_output = CompressedSubgraphs(input, key_vertex_id); 
            bool add = true;
            for (uint32_t set_vid : parents_) {
              std::vector<VertexID> new_set;
              uint32_t id = query_vertex_indices_[set_vid];
              intersect(*input.getSet(id), key_out_neighbors, &new_set);
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
          }
        }
      }
      inputs_idx_++;
      if (output_num >= batch_size) {
        break;
      }
    }
    return output_num;
  }
 
 protected:
  bool isInCandidates(VertexID key) {
    auto lb = std::lower_bound(candidates_->begin(), candidates_->end(), key);
    return lb != candidates_->end() && *lb == key;
  }

  const std::vector<CompressedSubgraphs>* inputs_;
  const Graph* data_graph_;
  uint32_t inputs_idx_;
};

}  // namespace circinus

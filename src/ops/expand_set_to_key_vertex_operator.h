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

#include <unordered_set>
#include <vector>

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"

namespace circinus {

class ExpandSetToKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandSetToKeyVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    while (input_index_ < current_inputs_->size()) {
      std::unordered_set<VertexID> visited;
      auto& input = (*current_inputs_)[input_index_];
      const auto& parent_match = input.getSet(query_vertex_indices_[parents_[0]]);
      for (VertexID vid : *parent_match) {
        const auto& out_neighbors = current_data_graph_->getOutNeighbors(vid);
        for (uint32_t i = 0; i < out_neighbors.second; ++i) {
          VertexID key_vertex_id = out_neighbors.first[i];
          if (visited.find(key_vertex_id) == visited.end()) {
            visited.insert(key_vertex_id);
            if (!isInCandidates(key_vertex_id)) {
              continue;
            }
            auto key_out_neighbors = current_data_graph_->getOutNeighbors(key_vertex_id);
            // TODO(by) hash key_out_neighbors
            CompressedSubgraphs new_output(input, key_vertex_id);
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
              outputs->emplace_back(std::move(new_output));
              output_num++;
              // TODO(by) bound by batch_size
            }
          }
        }
      }
      input_index_++;
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
};

}  // namespace circinus

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
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"

namespace circinus {

class FilterBase {
 protected:
  const QueryGraph* query_graph_;
  const Graph* data_graph_;

 public:
  FilterBase(const QueryGraph* query_graph, const Graph* data_graph)
      : query_graph_(query_graph), data_graph_(data_graph) {}

  void pruneByPivotVertices(QueryVertexID query_vertex, QueryVertexID pivot_vertex,
                            std::vector<std::vector<VertexID>>& candidates) {
    std::vector<QueryVertexID> pivot_vertices = {pivot_vertex};
    pruneByPivotVertices(query_vertex, pivot_vertices, candidates);
  }

  void pruneByPivotVertices(QueryVertexID query_vertex, const std::vector<QueryVertexID>& pivot_vertices,
                            std::vector<std::vector<VertexID>>& candidates) {
    LabelID query_vertex_label = query_graph_->getVertexLabel(query_vertex);
    uint32_t query_vertex_degree = query_graph_->getVertexOutDegree(query_vertex);
    uint32_t count = 0;
    std::unordered_map<VertexID, uint32_t> flag;

    for (QueryVertexID pivot_vertex : pivot_vertices) {
      for (VertexID v : candidates[pivot_vertex]) {
        const auto& v_nbrs = data_graph_->getOutNeighbors(v);
        for (uint32_t k = 0; k < v_nbrs.second; ++k) {
          VertexID v_nbr = v_nbrs.first[k];
          LabelID v_nbr_label = data_graph_->getVertexLabel(v_nbr);
          uint32_t v_nbr_degree = data_graph_->getVertexOutDegree(v_nbr);

          if (flag[v_nbr] == count && v_nbr_label == query_vertex_label && v_nbr_degree >= query_vertex_degree) {
            flag[v_nbr] += 1;
          }
        }
      }
      count += 1;
    }
    uint32_t new_candidate_size = 0;
    candidates[query_vertex].erase(std::remove_if(candidates[query_vertex].begin(), candidates[query_vertex].end(),
                                                  [&](const VertexID& u) { return flag[u] != count; }),
                                   candidates[query_vertex].end());
  }
};

}  // namespace circinus

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

#include <unordered_map>
#include <vector>

#include "ops/filters/ceci_filter.h"
#include "utils/utils.h"

namespace circinus {
CECIFilter::CECIFilter(const QueryGraph* query_graph, const Graph* data_graph, QueryVertexID start_vertex) {
  query_graph_ = query_graph;
  data_graph_ = data_graph;
  start_vertex_ = start_vertex;
  uint32_t query_vertices_num = query_graph->getNumVertices();
  bfs(query_graph, start_vertex, bfs_tree_, bfs_order_);
  std::vector<uint32_t> order_index(query_vertices_num);
  for (uint32_t i = 0; i < query_vertices_num; ++i) {
    QueryVertexID query_vertex = bfs_order_[i];
    order_index[query_vertex] = i;
  }

  for (uint32_t i = 0; i < query_vertices_num; ++i) {
    QueryVertexID u = bfs_order_[i];

    const auto& u_nbrs = query_graph->getOutNeighbors(u);
    for (uint32_t j = 0; j < u_nbrs.second; ++j) {
      QueryVertexID u_nbr = u_nbrs.first[j];
      if (u_nbr != bfs_tree_[u].parent_ && order_index[u_nbr] < order_index[u]) {
        bfs_tree_[u].bn_.emplace_back(u_nbr);
        bfs_tree_[u_nbr].fn_.emplace_back(u);
      }
    }
  }
}

}  // namespace circinus

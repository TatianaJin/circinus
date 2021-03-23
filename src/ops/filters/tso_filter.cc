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

#include "ops/filters/tso_filter.h"
#include "utils/utils.h"

namespace circinus {

TSOFilter::TSOFilter(const QueryGraph* query_graph, const Graph* data_graph, QueryVertexID start_vertex)
    : FilterBase(query_graph, data_graph), start_vertex_(start_vertex) {
  std::vector<QueryVertexID> bfs_order_;
  bfs(query_graph, start_vertex, tree_, bfs_order_);
  dfs(start_vertex, tree_, dfs_order_);
}

void TSOFilter::Filter(std::vector<std::vector<VertexID>>& candidates) {
  for (QueryVertexID query_vertex : dfs_order_) {
    if (query_vertex == start_vertex_) {
      continue;
    }
    TreeNode& node = tree_[query_vertex];
    pruneByPivotVertices(query_vertex, node.parent_, candidates);
  }

  for (auto it = dfs_order_.rbegin(); it != dfs_order_.rend(); ++it) {
    TreeNode& node = tree_[*it];
    pruneByPivotVertices(*it, node.children_, candidates);
  }
}

}  // namespace circinus

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

#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"

namespace circinus {

class CECIFilter {
 private:
  const QueryGraph* query_graph_;
  const Graph* data_graph_;
  QueryVertexID start_vertex_;
  std::vector<TreeNode> bfs_tree_;
  std::vector<QueryVertexID> bfs_order_;
  std::vector<QueryVertexID> level_offset_;
  uint32_t level_num_;

 public:
  CECIFilter(const QueryGraph* query_graph, const Graph* data_graph, QueryVertexID start_vertex);

  /** @returns The number of records that passed the filter and are added to output */
  void Filter(QueryVertexID query_vertex, std::vector<std::vector<VertexID>>& candidates);

 protected:
  void pruneByPivotVertices(QueryVertexID query_vertex, const std::vector<QueryVertexID>& pivot_vertices,
                            std::vector<std::vector<VertexID>>& candidates);
};

}  // namespace circinus

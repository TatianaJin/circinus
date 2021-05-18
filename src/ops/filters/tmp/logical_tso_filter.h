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

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters/logical_filter.h"

namespace circinus {

class LogicalTSOFilter : public LogicalFilter {
 public:
  LogicalTSOFilter(const QueryGraph* query_graph, const Graph* data_graph, const QueryVertexID query_vertex,
                   const QueryVertexID& pivot_vertex) {
    LogicalFilter(query_graph, data_graph, query_vertex, {pivot_vertex});
  }

  LogicalTSOFilter(const QueryGraph* query_graph, const Graph* data_graph, const QueryVertexID query_vertex,
                   const std::vector<QueryVertexID>& pivot_vertices) {
    LogicalFilter(query_graph, data_graph, query_vertex, pivot_vertices);
  }
};

}  // namespace circinus

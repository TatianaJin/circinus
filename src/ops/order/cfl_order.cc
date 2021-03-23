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

#include <queue>
#include <utility>
#include <vector>

#include "algorithms/k_core.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/order/cfl_order.h"

namespace circinus {

uint32_t CFLOrder::getStartVertex(const Graph* data_graph, const QueryGraph* query_graph,
                                  std::vector<uint32_t>& candidate_size) {
  double min_score = data_graph->getNumVertices();
  QueryVertexID start_vertex = 0;

  TwoCoreSolver two_core_solver;
  const auto& core_table = two_core_solver.get2CoreTable(query_graph);
  uint32_t core_size = two_core_solver.getCoreSize(core_table);
  for (QueryVertexID query_vertex = 0; query_vertex < query_graph->getNumVertices(); ++query_vertex) {
    if (core_size == 0 || core_table[query_vertex] > 1) {
      double degree = query_graph->getVertexOutDegree(query_vertex);
      double cur_score = candidate_size[query_vertex] / (double)degree;
      if (cur_score < min_score) {
        min_score = cur_score;
        start_vertex = query_vertex;
      }
    }
  }

  return start_vertex;
}

}  // namespace circinus

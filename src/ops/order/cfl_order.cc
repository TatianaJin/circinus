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
  auto rank_compare = [](std::pair<VertexID, double> l, std::pair<VertexID, double> r) { return l.second < r.second; };

  std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(rank_compare)>
      rank_queue(rank_compare);

  TwoCoreSolver two_core_solver;
  const auto& core_table = two_core_solver.get2CoreTable(query_graph);
  uint32_t core_size = two_core_solver.getCoreSize(core_table);
  for (QueryVertexID query_vertex = 0; query_vertex < query_graph->getNumVertices(); ++query_vertex) {
    if (core_size == 0 || core_table[query_vertex] > 1) {
      LabelID label = query_graph->getVertexLabel(query_vertex);
      double degree = query_graph->getVertexOutDegree(query_vertex);
      uint32_t frequency = data_graph->getNumVerticesByLabel(label);
      double rank = frequency / degree;
      rank_queue.push(std::make_pair(query_vertex, rank));
    }
  }

  // Keep the top-3.
  while (rank_queue.size() > 3) {
    rank_queue.pop();
  }

  QueryVertexID start_vertex = 0;
  double min_score = data_graph->getNumVertices() + 1;

  while (!rank_queue.empty()) {
    QueryVertexID query_vertex = rank_queue.top().first;
    double cur_score = candidate_size[query_vertex] / (double)query_graph->getVertexOutDegree(query_vertex);

    if (cur_score < min_score) {
      start_vertex = query_vertex;
      min_score = cur_score;
    }
    rank_queue.pop();
  }

  return start_vertex;
}

}  // namespace circinus

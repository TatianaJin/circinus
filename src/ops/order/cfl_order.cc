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

std::vector<QueryVertexID> CFLOrder::getTopThree(const GraphMetadata& metadata, const QueryGraph* q) {
  auto rank_compare = [](std::pair<QueryVertexID, double> l, std::pair<QueryVertexID, double> r) {
    return l.second < r.second;
  };
  std::priority_queue<std::pair<QueryVertexID, double>, std::vector<std::pair<QueryVertexID, double>>,
                      decltype(rank_compare)>
      rank_queue(rank_compare);
  TwoCoreSolver two_core_solver;
  const auto& core_table = two_core_solver.get2CoreTable(q);
  uint32_t core_size = two_core_solver.getCoreSize(core_table);

  for (QueryVertexID query_vertex = 0; query_vertex < q->getNumVertices(); ++query_vertex) {
    if (core_size == 0 || two_core_solver.isInCore(core_table, query_vertex)) {
      double degree = q->getVertexOutDegree(query_vertex);
      uint32_t frequency = metadata.getLabelFrequency(q->getVertexLabel(query_vertex));
      double rank = frequency / degree;
      rank_queue.emplace(query_vertex, rank);
    }
  }
  std::vector<QueryVertexID> ret(3);
  // Keep the top-3.
  while (rank_queue.size() > 3) {
    rank_queue.pop();
  }
  for (uint32_t i = 0; i < 3; ++i) {
    ret[i] = rank_queue.top().first;
    rank_queue.pop();
  }
  return ret;
}

QueryVertexID CFLOrder::getStartVertex(const std::vector<QueryVertexID>& query_vertices,
                                       const std::vector<VertexID>& candidate_size, const QueryGraph& query_graph,
                                       const GraphMetadata& metadata) {
  QueryVertexID start_vertex = 0;
  double min_score = metadata.getNumVertices() + 1;
  for (auto query_vertex : query_vertices) {
    double cur_score = candidate_size[query_vertex] / (double)query_graph.getVertexOutDegree(query_vertex);
    if (cur_score < min_score) {
      start_vertex = query_vertex;
      min_score = cur_score;
    }
  }
  return start_vertex;
}

QueryVertexID CFLOrder::getStartVertex(const Graph* data_graph, const QueryGraph* query_graph,
                                       const std::vector<uint32_t>& candidate_size) {
  GraphMetadata metadata(*data_graph);
  auto root_candidates = getTopThree(metadata, query_graph);
  return getStartVertex(root_candidates, {candidate_size.begin(), candidate_size.end()}, *query_graph, metadata);
}

}  // namespace circinus

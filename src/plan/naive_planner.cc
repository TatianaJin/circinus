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

#include "plan/naive_planner.h"

#include <queue>
#include <vector>

#include "algorithms/k_core.h"
#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/query_graph.h"
#include "plan/execution_plan.h"

namespace circinus {

bool NaivePlanner::hasValidCandidate() {
  for (auto& cardinality : *candidate_cardinality_) {
    if (cardinality < 1) {
      return false;
    }
  }
  return true;
}

ExecutionPlan* NaivePlanner::generatePlan() {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return nullptr;
  }
  vertex_cover_solver_.computeVertexCover();
  auto& covers = vertex_cover_solver_.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained

  // now we only consider a random MWVC from covers
  QueryVertexID v = 0;
  std::vector<QueryVertexID> cover;
  for (auto assignment : covers.front()) {
    if (assignment == 1) {
      cover.push_back(v);
    }
    ++v;
  }
  // start matching from the vertex with the smallest cardinality in cover
  auto start_vertex = selectStartingVertex(cover);

  TwoCoreSolver solver;
  auto& core_table = solver.get2CoreTable(query_graph_);
  auto order = generateMatchingOrder(query_graph_, core_table, start_vertex);
  plan_.populatePhysicalPlan(query_graph_, order, covers.front());
  return &plan_;
}

std::vector<QueryVertexID> NaivePlanner::generateMatchingOrder(const QueryGraph* g, const std::vector<int>& core_table,
                                                               QueryVertexID start_vertex) {
  std::vector<bool> visited(g->getNumVertices(), false);  // visited when pushed into the queue
  std::vector<QueryVertexID> order(g->getNumVertices());  // the order out of queue
  // if both in core or both not , prioritize the vertex with smaller cardinality; otherwise prioritize the core
  auto compare = [this, &core_table](QueryVertexID v1, QueryVertexID v2) {
    if (TwoCoreSolver::isInCore(core_table, v1)) {
      if (TwoCoreSolver::isInCore(core_table, v2)) {
        return candidate_cardinality_[v1] > candidate_cardinality_[v2];
      }
      return false;
    }
    if (TwoCoreSolver::isInCore(core_table, v2)) {
      return true;
    }
    return candidate_cardinality_[v1] > candidate_cardinality_[v2];
  };
  std::priority_queue<QueryVertexID, std::vector<QueryVertexID>, decltype(compare)> queue(compare);
  // initialization
  visited[start_vertex] = true;
  order.front() = start_vertex;
  uint32_t order_offset = 1;
  auto neighbors = g->getOutNeighbors(start_vertex);
  for (uint32_t i = 0; i < neighbors.second; ++i) {
    visited[neighbors.first[i]] = true;
    queue.push(neighbors.first[i]);
  }
  // repeatedly traverse to a best vertex
  while (!queue.empty()) {
    auto next = queue.top();
    queue.pop();
    order[order_offset++] = next;
    auto neighbors = g->getOutNeighbors(next);
    for (uint32_t i = 0; i < neighbors.second; ++i) {
      if (!visited[neighbors.first[i]]) {
        visited[neighbors.first[i]] = true;
        queue.push(neighbors.first[i]);
      }
    }
  }
  return order;
}

QueryVertexID NaivePlanner::selectStartingVertex(const std::vector<QueryVertexID>& cover) {
  // std::sort(cover.begin(), cover.end(), [this](QueryVertexID v1, QueryVertexID v2) { return
  // (*candidate_cardinality_)[v1] / query_graph_->getVertexOutDegree(v1) < (*candidate_cardinality_)[v2] /
  // query_graph_->getVertexOutDegree(v2); });

  QueryVertexID start_vertex;
  double min = (*candidate_cardinality_)[cover.front()];

  for (auto v : cover) {
    double score = (*candidate_cardinality_)[v] / query_graph_->getVertexOutDegree(v);
    if (score < min) {
      score = min;
      start_vertex = v;
    }
  }
  return start_vertex;
}

}  // namespace circinus

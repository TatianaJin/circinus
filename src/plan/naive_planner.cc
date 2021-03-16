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
#include <tuple>
#include <vector>

#include "algorithms/k_core.h"
#include "algorithms/minimum_weight_vertex_cover.h"
#include "algorithms/vertex_cover.h"
#include "graph/query_graph.h"
#include "plan/execution_plan.h"
#include "utils/hashmap.h"

namespace circinus {

bool NaivePlanner::hasValidCandidate() {
  for (auto& cardinality : *candidate_cardinality_) {
    if (cardinality < 1) {
      return false;
    }
  }
  return true;
}

ExecutionPlan* NaivePlanner::generatePlan(const std::vector<QueryVertexID>& use_order) {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return nullptr;
  }
  vertex_cover_solver_.computeVertexCover();
  auto& covers = vertex_cover_solver_.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto select_cover = covers[BnB::getSmallestCover(covers).first.front()];

  if (use_order.empty()) {
    TwoCoreSolver solver;
    auto& core_table = solver.get2CoreTable(query_graph_);
    // now we only consider a random smallest MWVC from covers
    QueryVertexID v = 0;
    std::vector<QueryVertexID> cover;
    for (auto assignment : select_cover) {
      if (assignment == 1) {
        cover.push_back(v);
      }
      ++v;
    }
    // start matching from the vertex with the smallest cardinality in cover
    auto start_vertex = selectStartingVertex(cover);

    matching_order_ = generateMatchingOrder(query_graph_, core_table, start_vertex);
    plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover);
  } else {
    matching_order_ = use_order;
    plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover);
  }
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
        return (*candidate_cardinality_)[v1] > (*candidate_cardinality_)[v2];
      }
      return false;
    }
    if (TwoCoreSolver::isInCore(core_table, v2)) {
      return true;
    }
    return (*candidate_cardinality_)[v1] > (*candidate_cardinality_)[v2];
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

  QueryVertexID start_vertex = cover.front();
  double min = (*candidate_cardinality_)[cover.front()];

  for (auto v : cover) {
    double score = (*candidate_cardinality_)[v] / query_graph_->getVertexOutDegree(v);
    DLOG(INFO) << "cover vertex " << v << " score " << score << " min " << min;
    if (score < min) {
      min = score;
      start_vertex = v;
    }
  }
  return start_vertex;
}

uint32_t NaivePlanner::analyzeDynamicCoreCoverEager(const std::vector<QueryVertexID>& use_order) {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return 0;
  }
  vertex_cover_solver_.computeVertexCover();
  auto& covers = vertex_cover_solver_.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto select_cover = covers[BnB::getSmallestCover(covers).first.front()];
  if (use_order.empty()) {  // generate a matching order if not specified
    TwoCoreSolver solver;
    auto& core_table = solver.get2CoreTable(query_graph_);
    // now we only consider a random smallest MWVC from covers
    QueryVertexID v = 0;
    std::vector<QueryVertexID> cover;
    for (auto assignment : select_cover) {
      if (assignment == 1) {
        cover.push_back(v);
      }
      ++v;
    }
    // start matching from the vertex with the smallest cardinality in cover
    auto start_vertex = selectStartingVertex(cover);
    matching_order_ = generateMatchingOrder(query_graph_, core_table, start_vertex);
  } else {
    matching_order_ = use_order;
  }

  uint32_t sum_delayed_steps = 0;
  auto cover = select_cover;
  unordered_set<QueryVertexID> existing_vertices;
  std::vector<uint32_t> vertex_order(matching_order_.size());
  existing_vertices.reserve(matching_order_.size());
  existing_vertices.insert(matching_order_.front());
  vertex_order[matching_order_.front()] = 0;
  // TODO(tatiana): use mwvc to decide the dynamic cover
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    auto v = matching_order_[i];
    vertex_order[v] = i;
    auto nb = query_graph_->getOutNeighbors(v);
    if (cover[v] == 1) {  // vertex v is in key, we eagerly check whether the current target can be put in key later
      bool delay_enumeration = true;
      // can be put in key later for a larger subquery if all the existing vertices connecting to v is in key
      for (uint32_t j = 0; j < nb.second; ++j) {
        if (existing_vertices.count(nb.first[j]) && cover[nb.first[j]] != 1) {
          delay_enumeration = false;
          break;
        }
      }
      cover[v] = !delay_enumeration;
    } else {
      for (uint32_t j = 0; j < nb.second; ++j) {  // restore the delayed vertex as key
        if (existing_vertices.count(nb.first[j]) && cover[nb.first[j]] != 1) {
          CHECK_EQ(select_cover[nb.first[j]], 1);
          cover[nb.first[j]] = 1;
          sum_delayed_steps += (i - vertex_order[nb.first[j]]);
          LOG(INFO) << nb.first[j] << " delayed from " << vertex_order[nb.first[j]] << " for "
                    << (i - vertex_order[nb.first[j]]) << " steps";
        }
      }
    }
    existing_vertices.insert(v);
  }
  return sum_delayed_steps;
}

std::tuple<uint32_t, uint32_t, uint32_t> NaivePlanner::analyzeDynamicCoreCoverMWVC() {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return std::tuple(0, 0, 0);
  }
  CHECK(!matching_order_.empty());

  // meters
  uint32_t sum_delayed_steps = 0;
  uint32_t sum_key_sizes = 0;
  uint32_t sum_key_to_set_changes = 0;

  auto subquery_vertices = matching_order_;
  auto candidate_cardinality = *candidate_cardinality_;

  // the static cover for the original query
  vertex_cover_solver_.computeVertexCover();
  auto& covers = vertex_cover_solver_.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto choice = BnB::getSmallestCover(covers);
  std::vector<bool> current_cover_assignment;
  current_cover_assignment.resize(matching_order_.size());
  std::stringstream ss;
  for (QueryVertexID i = 0; i < current_cover_assignment.size(); ++i) {
    current_cover_assignment[i] = (covers[choice.first.front()][i] == 1);
    ss << ' ' << (covers[choice.first.front()][matching_order_[i]] == 1);
  }
  LOG(INFO) << "cover" << ss.str() << " from " << choice.first.size();
  sum_key_sizes += choice.second;
  // compute the key sizes for static cover
  if (current_cover_assignment[matching_order_.front()]) {
    sum_delayed_steps += matching_order_.size() - 1;  // the first vertex appears from step 1
  }
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    auto v = matching_order_[i];
    if (current_cover_assignment[v]) {
      sum_delayed_steps += (matching_order_.size() - i);  // v is key for all the rest steps
    }
  }

  // use mwvc to decide the dynamic cover for each subquery i
  for (uint32_t i = matching_order_.size() - 2; i > 0; --i) {
    subquery_vertices.resize(i + 1);
    candidate_cardinality.resize(i + 1);
    // compute a vertex cover of subquery i
    auto subquery = query_graph_->getInducedSubgraph(subquery_vertices);
    WeightedBnB vc_solver(&subquery, candidate_cardinality);
    vc_solver.computeVertexCover();
    auto& covers = vc_solver.getBestCovers();
    CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
    auto choice = BnB::getSmallestCover(covers);
    auto select_cover = covers[choice.first.front()];
    sum_key_sizes += choice.second;
    // for existing vertices 0 ~ i, check the change in key/set assignment
    std::stringstream ss;
    for (QueryVertexID j = 0; j <= i; ++j) {
      if (select_cover[j] == 1) {  // now assign j to key
        // add to sum_key_to_set_changes if change from key to set
        sum_key_to_set_changes += (!current_cover_assignment[matching_order_[j]]);
      }
      ss << ' ' << (select_cover[j] == 1);
      current_cover_assignment[matching_order_[j]] = (select_cover[j] == 1);
    }
    LOG(INFO) << "cover" << ss.str() << " from " << choice.first.size();
  }

  sum_delayed_steps -= sum_key_sizes;
  return std::tuple(sum_key_sizes, sum_delayed_steps, sum_key_to_set_changes);
}

}  // namespace circinus

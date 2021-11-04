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

#include "algorithms/vertex_cover.h"

#include <cassert>
#include <deque>
#include <vector>

#include "graph/query_graph.h"

namespace circinus {

int BnB::countAssignment(const std::vector<int>& assignment) {
  int count = 0;
  for (uint32_t i = 0; i < assignment.size(); ++i) {
    count += (assignment[i] == 1);
  }
  return count;
}

std::deque<QueryEdge> BnB::getEdgeList(const QueryGraph& g) {
  std::deque<QueryEdge> list;
  for (QueryVertexID i = 0; i < g.getNumVertices(); ++i) {
    auto neighbors = g.getOutNeighbors(i);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (i < neighbors.first[j]) {
        list.emplace_back(i, neighbors.first[j]);
      }
    }
  }
  return list;
}

/** Use a 2-approximate algorithm to calculate a lower bound for a vertex cover of the uncovered edges. */
double BnB::matchingLB(std::deque<QueryEdge> uncovered_edges) {
  int count = 0;
  while (!uncovered_edges.empty()) {
    ++count;  // 1/2 of the selected vertices, each time we add one to count as we select two vertices
    int src = uncovered_edges.front().src;
    int dst = uncovered_edges.front().dst;
    for (int i = uncovered_edges.size() - 1; i >= 0; --i) {
      if (uncovered_edges[i].hasEndVertex(src) || uncovered_edges[i].hasEndVertex(dst)) {
        uncovered_edges.erase(uncovered_edges.begin() + i);
      }
    }
  }
  return count;
}

uint32_t BnB::countUncoveredNeighbors(int v, const std::deque<QueryEdge>& uncovered_edges) {
  uint32_t count = 0;
  for (uint32_t i = 0; i < uncovered_edges.size(); ++i) {
    count += uncovered_edges[i].hasEndVertex(v);
  }
  return count;
}

void BnB::dfs(std::vector<int> assignment, std::deque<QueryEdge> uncovered_edges, int depth) {
  // if the time limit is reached, stop searching and return the result
  if (((double)(clock() - start_time_) / CLOCKS_PER_SEC) > cutoff_time_) {
    return;
  }

  auto current_count = countAssignment(assignment);
  // prune when the current count already has more vertices than the best solution found so far
  if (current_count > best_cover_size_) {
    return;
  }

  // prune when the current count + lower bound to cover remaining edges > best result so far
  if (current_count + matchingLB(uncovered_edges) > best_cover_size_) {
    return;
  }

  // check infeasible assignment: at least one edge with both ends assigned to be not in the vertex cover
  for (uint32_t i = 0; i < uncovered_edges.size(); ++i) {
    if (assignment[uncovered_edges[i].src] == 0 && assignment[uncovered_edges[i].dst] == 0) {
      return;
    }
  }

  // when a vertex cover is found
  if (uncovered_edges.empty()) {
    if (current_count < best_cover_size_) {  // a cover better than current best covers is found
      best_cover_size_ = current_count;
      history_cover_size_.push_back(current_count);
      history_time_.push_back(((double)(clock() - start_time_) / CLOCKS_PER_SEC));
      best_covers_.resize(1);
      best_covers_.front().swap(assignment);
    } else if (current_count == best_cover_size_) {  // store all best covers
      best_covers_.emplace_back(std::move(assignment));
    }
    return;
  }

  // find the branch vertex that has most uncovered edges
  uint32_t branching_vertex = 0;
  uint32_t max_uncovered_neighbors = 0;
  for (uint32_t i = 0; i < graph_->getNumVertices(); ++i) {
    if (assignment[i] == -1) {
      auto tmp = countUncoveredNeighbors(i, uncovered_edges);
      if (tmp > max_uncovered_neighbors) {
        max_uncovered_neighbors = tmp;
        branching_vertex = i;
      }
    }
  }

  // first branch: not selecting branching_vertex
  assignment[branching_vertex] = 0;
  dfs(assignment, uncovered_edges, depth + 1);

  // second branch: selecting branching_vertex
  for (int i = uncovered_edges.size() - 1; i >= 0; --i) {
    if (uncovered_edges[i].hasEndVertex(branching_vertex)) uncovered_edges.erase(uncovered_edges.begin() + i);
  }
  assignment[branching_vertex] = 1;
  dfs(assignment, uncovered_edges, depth + 1);
}

}  // namespace circinus

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

#include "algorithms/minimum_weight_vertex_cover.h"

#include <cassert>
#include <deque>
#include <utility>
#include <vector>

#include "graph/query_graph.h"

namespace circinus {

double WeightedBnB::computeAssignmentWeight(const std::vector<int>& assignment,
                                            const std::vector<double>& vertex_weights) {
  double weight = 0;
  for (uint32_t i = 0; i < assignment.size(); ++i) {
    if (assignment[i] == 1) {
      weight += vertex_weights[i];
    }
  }
  return weight;
}

std::deque<QueryEdge> WeightedBnB::getEdgeList(const QueryGraph& g) {
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

// FIXME(tatiana): check the correctness of using this lower-bound
/** Use a 2-approximate algorithm to calculate a lower bound for a vertex cover of the uncovered edges. */
double WeightedBnB::matchingLB(std::deque<QueryEdge> uncovered_edges, std::vector<double> untight_weights) {
  double weight = 0;
  while (!uncovered_edges.empty()) {
    // ++count;  // 1/2 of the selected vertices, each time we add one to count as we select two vertices
    int src = uncovered_edges.front().src;
    int dst = uncovered_edges.front().dst;
    if (untight_weights[src] < untight_weights[dst]) {
      // make src vertex tight by assigning this edge a price of untight_weights[src]
      weight += vertex_weights_[src];
      untight_weights[dst] -= untight_weights[src];
      untight_weights[src] = 0;
      for (int i = uncovered_edges.size() - 1; i >= 0; --i) {
        if (uncovered_edges[i].hasEndVertex(src)) {
          uncovered_edges.erase(uncovered_edges.begin() + i);
        }
      }
    } else if (untight_weights[src] > untight_weights[dst]) {
      // make dst vertex tight by assigning this edge a price of untight_weights[dst]
      weight += vertex_weights_[dst];
      untight_weights[src] -= untight_weights[dst];
      untight_weights[dst] = 0;
      for (int i = uncovered_edges.size() - 1; i >= 0; --i) {
        if (uncovered_edges[i].hasEndVertex(dst)) {
          uncovered_edges.erase(uncovered_edges.begin() + i);
        }
      }
    } else {
      // make both ends tight
      weight += vertex_weights_[src] + vertex_weights_[dst];
      untight_weights[src] = 0;
      untight_weights[dst] = 0;
      for (int i = uncovered_edges.size() - 1; i >= 0; --i) {
        if (uncovered_edges[i].hasEndVertex(src) || uncovered_edges[i].hasEndVertex(dst)) {
          uncovered_edges.erase(uncovered_edges.begin() + i);
        }
      }
    }
  }
  return weight / 2;
}

uint32_t WeightedBnB::countUncoveredNeighbors(int v, const std::deque<QueryEdge>& uncovered_edges) {
  uint32_t count = 0;
  for (uint32_t i = 0; i < uncovered_edges.size(); ++i) {
    count += uncovered_edges[i].hasEndVertex(v);
  }
  return count;
}

void WeightedBnB::dfs(std::vector<int> assignment, std::deque<QueryEdge> uncovered_edges, int depth) {
  // if the time limit is reached, stop searching and return the result
  if (((double)(clock() - start_time_) / CLOCKS_PER_SEC) > cutoff_time_) {
    return;
  }

  auto current_weight = computeAssignmentWeight(assignment, vertex_weights_);
  // prune when the current count already has more vertices than the best solution found so far
  if (current_weight > best_cover_weight_) {
    return;
  }

  // prune when the current count + lower bound to cover remaining edges > best result so far
  if (current_weight + matchingLB(uncovered_edges, vertex_weights_) > best_cover_weight_) {
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
    if (best_covers_.size() < best_buffer_size_) {
      best_covers_.push_back(std::move(assignment));
      best_cover_weights_.push_back(current_weight);
    } else {  // a cover better than current best covers is found
      uint32_t idx = 0;
      for (uint32_t i = 1; i < best_buffer_size_; ++i) {
        if (best_cover_weights_[i] > best_cover_weights_[idx]) {
          idx = i;
        }
      }
      if (best_cover_weights_[idx] > current_weight) {
        best_cover_weights_[idx] = current_weight;
        history_cover_size_.push_back(current_weight);
        history_time_.push_back(((double)(clock() - start_time_) / CLOCKS_PER_SEC));
        best_covers_[idx].swap(assignment);
      }
    }
    return;
  }

  // find the branch vertex that has most uncovered edges
  uint32_t branching_vertex = 0;
  uint32_t max_uncovered_neighbors = 0;
  for (uint32_t i = 0; i < graph_->getNumVertices(); ++i) {
    if (assignment[i] == -1) {
      double tmp = countUncoveredNeighbors(i, uncovered_edges);
      tmp /= (vertex_weights_[i] + 1);  // add one to the denominator to avoid divide by zero
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
    if (uncovered_edges[i].hasEndVertex(branching_vertex)) {
      uncovered_edges.erase(uncovered_edges.begin() + i);
    }
  }
  assignment[branching_vertex] = 1;
  dfs(assignment, uncovered_edges, depth + 1);
}

}  // namespace circinus

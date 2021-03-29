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

#include <cassert>
#include <cstdint>
#include <ctime>
#include <deque>
#include <vector>

#include "graph/query_graph.h"

namespace circinus {

/**
 * A branch-and-bound algorithm to find exact minimum vertex cover.
 */
class WeightedBnB {
 private:
  const QueryGraph* graph_;
  const std::vector<double>& vertex_weights_;
  clock_t start_time_;
  double cutoff_time_;                         // time limit for the algorithm
  std::vector<std::vector<int>> best_covers_;  // the best vertex covers found so far
  double best_cover_weight_;                   // the size of the best vertex cover found so far
  double elapsed_time_;

  // for analyzing the algorithm efficiency
  // TODO(tatiana): put for debug mode only?
  std::vector<int> history_cover_size_;
  std::vector<double> history_time_;

 public:
  /**
   * @param g The query graph for which we find the minimum vertex cover(s).
   * @param cutoff_time The time limit of the execution. The best cover found within cutoff_time will be returned
   */
  WeightedBnB(const QueryGraph* g, const std::vector<double>& vertex_weights, double cutoff_time = 10)
      : graph_(g), vertex_weights_(vertex_weights), cutoff_time_(cutoff_time) {
    best_cover_weight_ = 0;
    for (auto w : vertex_weights_) {
      best_cover_weight_ += w;
    }
  }

  /** Compute minimum vertex cover(s) within cutoff_time_ */
  inline void computeVertexCover() {
    if (!best_covers_.empty()) return;
    start_time_ = clock();
    dfs(std::vector<int>(graph_->getNumVertices(), -1), getEdgeList(*graph_), 0);
    elapsed_time_ = ((double)(clock() - start_time_)) / CLOCKS_PER_SEC;
  }

  // getters
  inline const auto& getBestCovers() const { return best_covers_; }
  inline double getBestObjective() const { return best_cover_weight_; }
  /** * @returns Time used to find the first best cover */
  inline double getTimeToBest() const { return history_time_.back(); }
  /** * @returns Time used to find the best covers */
  inline double getElapsedTime() const { return elapsed_time_; }

 private:
  static std::deque<QueryEdge> getEdgeList(const QueryGraph& g);
  static double computeAssignmentWeight(const std::vector<int>& assignment, const std::vector<double>& vertex_weights);
  /** Use a 2-approximate algorithm to calculate a lower bound for a vertex cover of the uncovered edges. */
  static uint32_t countUncoveredNeighbors(int v, const std::deque<QueryEdge>& uncovered_edges);

  double matchingLB(std::deque<QueryEdge> uncovered_edges, std::vector<double> untight_weights);

  // Branch and bound
  // assignment -1: undetermined, 0: not in vertex cover, 1: in vertex cover
  void dfs(std::vector<int> assignment, std::deque<QueryEdge> uncovered_edges, int depth);
};

}  // namespace circinus

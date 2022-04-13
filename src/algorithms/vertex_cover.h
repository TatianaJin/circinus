#pragma once

#include <cassert>
#include <cstdint>
#include <ctime>
#include <deque>
#include <utility>
#include <vector>

#include "graph/query_graph.h"

namespace circinus {

/**
 * A branch-and-bound algorithm to find exact minimum vertex cover.
 */
class BnB {
 private:
  const QueryGraph* graph_;
  clock_t start_time_;
  double cutoff_time_;                         // time limit for the algorithm
  std::vector<std::vector<int>> best_covers_;  // the best vertex covers found so far
  int best_cover_size_;                        // the size of the best vertex cover found so far
  double elapsed_time_;

  // for analyzing the algorithm efficiency
  std::vector<int> history_cover_size_;
  std::vector<double> history_time_;

 public:
  /**
   * @param g The query graph for which we find the minimum vertex cover(s).
   * @param cutoff_time The time limit of the execution. The best cover found within cutoff_time will be returned
   */
  explicit BnB(const QueryGraph* g, double cutoff_time = 5)
      : graph_(g), cutoff_time_(cutoff_time), best_cover_size_(g->getNumVertices()) {}

  /** Compute minimum vertex cover(s) within cutoff_time_ */
  inline void computeVertexCover() {
    start_time_ = clock();
    dfs(std::vector<int>(graph_->getNumVertices(), -1), getEdgeList(*graph_), 0);
    elapsed_time_ = ((double)(clock() - start_time_)) / CLOCKS_PER_SEC;
  }

  // getters
  inline const auto& getBestCovers() const { return best_covers_; }
  inline int getBestCoverSize() const { return best_cover_size_; }
  /** * @returns Time used to find the first best cover */
  inline double getTimeToBest() const { return history_time_.back(); }
  /** * @returns Time used to find the best covers */
  inline double getElapsedTime() const { return elapsed_time_; }

  /** @returns The index of the smallest cover in cover_assignments and the cover size */
  static inline std::pair<std::vector<uint32_t>, uint32_t> getSmallestCover(
      const std::vector<std::vector<int>>& cover_assignments) {
    DCHECK(!cover_assignments.empty());
    uint32_t min = cover_assignments.front().size();
    std::vector<uint32_t> ret;
    for (uint32_t i = 0; i < cover_assignments.size(); ++i) {
      uint32_t size = countAssignment(cover_assignments[i]);
      if (size < min) {
        min = size;
        ret.resize(1);
        ret.front() = i;
      } else if (size == min) {
        ret.push_back(i);
      }
    }
    return std::make_pair(ret, min);
  }

 private:
  static std::deque<QueryEdge> getEdgeList(const QueryGraph& g);
  static int countAssignment(const std::vector<int>& assignment);
  /** Use a 2-approximate algorithm to calculate a lower bound for a vertex cover of the uncovered edges. */
  static double matchingLB(std::deque<QueryEdge> uncovered_edges);
  static uint32_t countUncoveredNeighbors(int v, const std::deque<QueryEdge>& uncovered_edges);

  // Branch and bound
  // assignment -1: undetermined, 0: not in vertex cover, 1: in vertex cover
  void dfs(std::vector<int> assignment, std::deque<QueryEdge> uncovered_edges, int depth);
};

}  // namespace circinus

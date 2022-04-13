#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "algorithms/intersect.h"
#include "graph/query_graph.h"
#include "graph/query_neighbor_set.h"
#include "graph/types.h"

namespace circinus {

template <uint32_t K, bool count_label = false>
class KCycle {
 private:
  const QueryGraph* graph_;
  std::array<uint32_t, K - 2> cycle_counts_;  // k-3: count of k-cycle (induced subgraph)
  unordered_map<std::string, uint32_t> labeled_cycle_counts_;

 public:
  explicit KCycle(const QueryGraph* graph) : graph_(graph) {
    static_assert(K > 3);
    cycle_counts_.fill(0);
    auto n_qvs = graph_->getNumVertices();
    std::array<QueryVertexID, K> cycle;
    std::vector<bool> cycle_mask(graph_->getNumVertices(), false);  // neighbors of head vertex
    std::vector<bool> filter(graph_->getNumVertices(), false);      // union of neighbors of middle vertices
    for (cycle[0] = 0; cycle[0] < n_qvs; ++cycle[0]) {
      auto nbrs = graph_->getOutNeighbors(cycle[0]);
      for (auto end = nbrs.first + nbrs.second - 1; end >= nbrs.first; --end) {
        if (*end <= cycle[0]) break;  // consider only neighbors with larger id
        cycle_mask[*end] = true;
      }
      for (auto end = nbrs.first + nbrs.second - 1; end >= nbrs.first; --end) {
        if (*end <= cycle[0]) break;  // consider only neighbors with larger id
        cycle[1] = *end;
        dfs<2>(cycle, cycle_mask, filter);
      }
      cycle_mask.assign(graph_->getNumVertices(), false);
    }
  }

  uint32_t getCycleCount(uint32_t k) const {
    CHECK_GE(k, 3);
    DCHECK_LE(k, K);
    return cycle_counts_[k - 3];
  }

  const auto& getLabeledCycleCounts() const { return labeled_cycle_counts_; }

 private:
  /** @param cycle A simple path that is a induced subgraph which may extend into a cycle. */
  template <uint32_t level>
  inline void dfs(std::array<QueryVertexID, K>& cycle, const std::vector<bool>& cycle_mask,
                  const std::vector<bool>& filter) {
    auto tail_vertex = cycle[level - 1];
    auto nbrs = graph_->getOutNeighbors(tail_vertex);
    std::vector<bool> new_filter;
    if
      constexpr(level + 1 != K) {  // build new filter by treating the current tail as a middle vertex
        new_filter = filter;
        for (uint32_t nb = 0; nb < nbrs.second; ++nb) {
          new_filter[nbrs.first[nb]] = true;
        }
      }

    // the new vertex should be a neighbor of cycle[level - 1] but not a neighbor of any vertex in the middle
    for (auto end = nbrs.first + nbrs.second - 1; end >= nbrs.first; --end) {
      if (*end <= tail_vertex) break;  // consider only neighbors with larger id
      if (filter[*end]) continue;
      // two cases: find a simple cycle, or find a simple path
      if (cycle_mask[*end]) {  // if a cycle is found
        ++cycle_counts_[level - 2];
        if (level > 3 && count_label) {
          std::vector<LabelID> labels(level);
          for (uint32_t i = 0; i < level; ++i) {
            labels[i] = graph_->getVertexLabel(cycle[i]);
          }
          std::sort(labels.begin(), labels.end());
          std::stringstream ss;
          ss << "cycle";
          for (auto l : labels) {
            ss << '|' << l;
          }
          ++labeled_cycle_counts_[ss.str()];
        }
      } else {
        if
          constexpr(level + 1 == K) {
            return;  // stop recursion
          }
        cycle[level] = *end;                                         // new tail
        dfs<std::min(level + 1, K)>(cycle, cycle_mask, new_filter);  // stop template recursion at K
      }
    }
  }
};  // class KCycle

}  // namespace circinus

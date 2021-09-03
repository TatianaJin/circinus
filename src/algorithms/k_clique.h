// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
// the License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

#pragma once

#include <vector>

#include "algorithms/intersect.h"
#include "graph/query_graph.h"
#include "graph/query_neighbor_set.h"
#include "graph/types.h"

namespace circinus {

template <uint32_t K, bool count_label = false>
class KClique {
 private:
  const QueryGraph* graph_;
  std::array<uint32_t, K - 2> clique_counts_;  // k-3: count of k-clique
  unordered_map<std::string, uint32_t> labeled_clique_counts_;

 public:
  explicit KClique(const QueryGraph* graph) : graph_(graph) {
    static_assert(K > 2);
    clique_counts_.fill(0);
    auto n_qvs = graph_->getNumVertices();
    std::array<QueryVertexID, K> clique;
    for (clique[0] = 0; clique[0] < n_qvs; ++clique[0]) {
      auto nbrs = graph_->getOutNeighbors(clique[0]);
      for (auto end = nbrs.first + nbrs.second - 1; end >= nbrs.first; --end) {
        if (*end <= clique[0]) break;  // consider only neighbors with larger id
        clique[1] = *end;
        dfs<2>(clique);
      }
    }
  }

  uint32_t getCliqueCount(uint32_t k) const {
    if (k == 2) return graph_->getNumEdges();
    DCHECK_LE(k, K);
    return clique_counts_[k - 3];
  }

  const auto& getLabeledCliqueCounts() const { return labeled_clique_counts_; }

 private:
  template <uint32_t level>
  inline void dfs(std::array<QueryVertexID, K>& clique) {
    // count for each found clique
    if
      constexpr(level > 2) {
        ++clique_counts_[level - 3];
        if
          constexpr(count_label) {
            std::vector<LabelID> labels(level);
            for (uint32_t i = 0; i < level; ++i) {
              labels[i] = graph_->getVertexLabel(clique[i]);
            }
            std::sort(labels.begin(), labels.end());
            std::stringstream ss;
            ss << "clique";
            for (auto l : labels) {
              ss << '|' << l;
            }
            ++labeled_clique_counts_[ss.str()];
          }
      }
    if
      constexpr(level == K) return;

    std::array<QueryNeighborSet, level> neighbor_sets;
    // the new vertex should be common neighbors with id larger than that of clique[level - 1]
    for (uint32_t i = 0; i < level; ++i) {
      auto nbrs = graph_->getOutNeighbors(clique[i]);
      DCHECK_NE(nbrs.second, 0);  // connected graph min degree >= 1
      uint32_t nb = 0;
      while (nb < nbrs.second && nbrs.first[nb] < clique[level - 1]) {
        ++nb;
      }
      if (nb == nbrs.second) {
        return;  // no extension
      }
      neighbor_sets[i].addRange(nbrs.first + nb, nbrs.first + nbrs.second);
    }
    std::sort(neighbor_sets.begin(), neighbor_sets.end());

    // intersection for common neighbors
    std::vector<QueryVertexID> extensions;
    intersect(neighbor_sets[0], neighbor_sets[1], &extensions);
    for (uint32_t i = 2; i < level; ++i) {
      intersectInplace(extensions, neighbor_sets[i], &extensions);
    }
    // recursive search
    for (auto extension : extensions) {
      clique[level] = extension;
      dfs<std::min(level + 1, K)>(clique);  // stop template recursion at K
    }
  }
};  // class KClique

}  // namespace circinus

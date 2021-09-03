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

/** A k-star is a 1-to-k bipartite graph.
 * @tparam K: the value of k for k-star.
 */
template <uint32_t K>
class KStar {
 private:
  const QueryGraph* graph_;
  // k-2: count of k-star (induced subgraph, a k-star is not counted if part of a (k+1)-star)
  std::array<uint32_t, K - 1> star_counts_;

 public:
  explicit KStar(const QueryGraph* graph) : graph_(graph) {
    static_assert(K > 1);
    using VID = QueryVertexID;
    star_counts_.fill(0);

    auto n_qvs = graph_->getNumVertices();
    for (VID v = 0; v < n_qvs; ++v) {
      auto nbrs = graph_->getOutNeighbors(v);
      uint32_t k = 0;
      // TODO(tatiana): now only consider a center node with degree-1 neighbors, later may consider CFL-like
      // decomposition to count forest paths around a center node
      for (uint32_t nb = 0; nb < nbrs.second; ++nb) {
        if (graph_->getVertexOutDegree(nbrs.first[nb]) == 1) {
          ++k;
        }
      }
      if (k >= 2) {
        ++star_counts_[k - 2];
      }
    }
  }

  uint32_t getStarCount(uint32_t k) const {
    CHECK_GE(k, 2);
    DCHECK_LE(k, K);
    return star_counts_[k - 2];
  }
};  // class KStar

}  // namespace circinus

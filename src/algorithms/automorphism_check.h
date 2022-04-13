#pragma once

#include <map>
#include <set>
#include <utility>
#include <vector>

#include "bliss/graph.hh"

#include "algorithms/partial_order.h"
#include "graph/query_graph.h"
#include "graph/types.h"

namespace circinus {

class AutomorphismCheck {
 private:
  bliss::Graph bliss_q_;
  const QueryGraph* const q_;
  std::vector<QueryVertexID> mapping_;   // original to sorted qv
  std::vector<QueryVertexID> sorted_v_;  // sorted to original qv

 public:
  explicit AutomorphismCheck(QueryGraph& q);

  /** Each returned pair means a partial order first < second */

  PartialOrder getPartialOrder();

 private:
  static inline bool compareNeighbors(const std::pair<const QueryVertexID*, uint32_t>& nbrs1,
                                      const std::pair<const QueryVertexID*, uint32_t>& nbrs2) {
    DCHECK_EQ(nbrs1.second, nbrs2.second);
    for (uint32_t i = 0; i < nbrs1.second; ++i) {
      if (nbrs1.first[i] != nbrs2.first[i]) {
        return nbrs1.first[i] > nbrs2.first[i];
      }
    }
    return true;
  }

  PartialOrder orderVertices(const std::vector<std::pair<QueryVertexID, QueryVertexID>>& sorted_conds);

  /** @returns Permutations of vertices */
  inline std::vector<std::vector<uint32_t>> getAutomorphisms() {
    std::vector<std::vector<uint32_t>> result;
    bliss::Stats s;
    bliss_q_.find_automorphisms(s,
                                [](void* resultp, uint32_t n_vertices, const uint32_t* aut) {
                                  auto res = ((std::vector<std::vector<uint32_t>>*)resultp);
                                  res->push_back(std::vector<uint32_t>(aut, aut + n_vertices));
                                },
                                &result);
    return result;
  }

  std::map<uint32_t, std::set<uint32_t>> getAEquivalenceClasses(const std::vector<std::vector<uint32_t>>& A);
};

}  // namespace circinus

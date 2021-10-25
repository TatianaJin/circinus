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

#include "algorithms/vertex_equivalence.h"

#include <unordered_set>
#include <utility>

#include "algorithms/partial_order.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "utils/flags.h"
#include "utils/hashmap.h"

namespace circinus {

VertexEquivalence::VertexEquivalence(const QueryGraph& q, PartialOrder* po) : q_(&q), po_(po) {
  auto n_qvs = q.getNumVertices();
  for (QueryVertexID u1 = 0; u1 < n_qvs; ++u1) {  // find all equivalent to u1
    auto u1_nbrs = q.getOutNeighbors(u1);
    for (QueryVertexID u2 = u1 + 1; u2 < n_qvs; ++u2) {
      // check label equivalence
      if (q.getVertexLabel(u1) != q.getVertexLabel(u2)) continue;
      auto u2_nbrs = q.getOutNeighbors(u2);
      if (u1_nbrs.second == u2_nbrs.second) {
        bool equivalent = true;
        // strictly equal neighborhood, do not consider two adjacent vertices now
        for (uint32_t i = 0; i < u1_nbrs.second; ++i) {
          if (u1_nbrs.first[i] != u2_nbrs.first[i]) {
            equivalent = false;
            break;
          }
        }
        if (!equivalent) continue;
        // the partial order constraints need also be equal, except for the constraint between u1 and u2
        if (po != nullptr && !po->constraints.empty()) {
          auto u1_cons = po->constraints.find(u1);
          auto u2_cons = po->constraints.find(u2);
          auto end = po->constraints.end();
          if (u1_cons == end) {
            equivalent = (u2_cons == end);
          } else if (u2_cons == end) {
            equivalent = false;
          } else {
            equivalent = setsEqualExcept(u1_cons->second.first, u2_cons->second.first, u1, u2) &&
                         setsEqualExcept(u1_cons->second.second, u2_cons->second.second, u1, u2);
          }
        }
        if (equivalent) {
          equivalent_pairs_.emplace(u1, u2);
        }
      }
    }
  }
  if (verbosePlannerLog()) {
    std::stringstream ss;
    for (auto& p : equivalent_pairs_) {
      ss << ' ' << p.first << '=' << p.second;
    }
    LOG(INFO) << "Equivalent pairs" << ss.str();
  }
}

bool VertexEquivalence::setsEqualExcept(const unordered_set<QueryVertexID>& set1,
                                        const unordered_set<QueryVertexID>& set2, QueryVertexID u1, QueryVertexID u2) {
  if (set1.size() == set2.size()) {
    for (auto v : set1) {
      if (set2.count(v) == 0) return false;
    }
    return true;
  }
  if (set1.size() < set2.size()) {
    if (set1.size() + 1 != set2.size()) return false;
    if (set2.count(u1) == 0) return false;
    for (auto v : set1) {
      if (set2.count(v) == 0) return false;
    }
    return true;
  }
  if (set2.size() + 1 != set1.size()) return false;
  if (set1.count(u2) == 0) return false;
  for (auto v : set2) {
    if (set1.count(v) == 0) return false;
  }
  return true;
}

}  // namespace circinus

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

#include "algorithms/automorphism_check.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

#include "bliss/graph.hh"

#include "graph/query_graph.h"
#include "graph/types.h"
#include "utils/flags.h"

namespace circinus {

AutomorphismCheck::AutomorphismCheck(QueryGraph& q)
    : q_(&q), mapping_(q_->getNumVertices()), sorted_v_(q_->getNumVertices()) {
  // construct bliss graph with vertices sorted by degree
  QueryVertexID n_qvs = q.getNumVertices();

  std::iota(sorted_v_.begin(), sorted_v_.end(), 0);
  std::sort(sorted_v_.begin(), sorted_v_.end(), [&q, this](QueryVertexID a, QueryVertexID b) {
    return q.getVertexOutDegree(a) > q.getVertexOutDegree(b) ||
           (q.getVertexOutDegree(a) == q.getVertexOutDegree(b) &&
            compareNeighbors(q.getOutNeighbors(a), q.getOutNeighbors(b)));
  });
  for (uint32_t i = 0; i < n_qvs; ++i) {
    mapping_[sorted_v_[i]] = i;
    bliss_q_.add_vertex(q_->getVertexLabel(sorted_v_[i]));
  }
  for (uint32_t i = 0; i < n_qvs; ++i) {
    auto nbrs = q.getOutNeighbors(i);
    LOG(INFO) << "vertex " << i << " mapping " << mapping_[i] << " first neighbor " << nbrs.first[0];
    for (uint32_t nbr_i = 0; nbr_i < nbrs.second && i > nbrs.first[nbr_i]; ++nbr_i) {
      DLOG(INFO) << "edge " << mapping_[i] << ' ' << mapping_[nbrs.first[nbr_i]] << " from " << i << " "
                 << nbrs.first[nbr_i];
      bliss_q_.add_edge(mapping_[i], mapping_[nbrs.first[nbr_i]]);
    }
  }
}

std::vector<std::pair<QueryVertexID, QueryVertexID>> AutomorphismCheck::getPartialOrder() {
  std::vector<std::pair<QueryVertexID, QueryVertexID>> result;
  auto automorphisms = getAutomorphisms();
  if (verbosePlannerLog()) {
    std::stringstream ss;
    for (auto& permutation : automorphisms) {
      ss << std::endl;
      for (auto v : permutation) ss << " " << v;
    }
    LOG(INFO) << "automorphisms:" << ss.str();
  }

  auto eclasses = getAEquivalenceClasses(automorphisms);
  std::vector<std::pair<QueryVertexID, QueryVertexID>> conds;
  // Don't pick largest: want conditions on what will likely be the core vertices first, so that there are no
  // conditions across sibling groups.
  auto eclass_it = std::find_if(eclasses.cbegin(), eclasses.cend(), [](auto&& a1) { return a1.second.size() > 1; });
  while (eclass_it != eclasses.cend() && eclass_it->second.size() > 1) {
    const auto& eclass = eclass_it->second;

    if (verbosePlannerLog()) {
      std::stringstream ss;
      for (auto&& s : eclasses) {
        ss << std::endl << s.first << ": ";
        for (QueryVertexID v : s.second) ss << v << " ";
      }
      LOG(INFO) << "equivalence classes:" << ss.str();
    }

    QueryVertexID n0 = *eclass.cbegin();
    for (auto&& f : automorphisms) {
      QueryVertexID min = *std::min_element(std::next(eclass.cbegin()), eclass.cend(),
                                            [&f](QueryVertexID n, QueryVertexID m) { return f[n] < f[m]; });
      conds.emplace_back(n0, min);
      if (verbosePlannerLog()) {
        LOG(INFO) << "add partial order condition " << n0 << '<' << min;
      }
    }

    automorphisms.erase(std::remove_if(automorphisms.begin(), automorphisms.end(),
                                       [n0](auto& permutation) { return permutation[n0] != n0; }),
                        automorphisms.end());
    if (verbosePlannerLog()) {
      std::stringstream ss;
      for (auto& permutation : automorphisms) {
        ss << std::endl;
        for (auto v : permutation) ss << " " << v;
      }
      LOG(INFO) << "automorphisms:" << ss.str();
    }

    eclasses = getAEquivalenceClasses(automorphisms);
    eclass_it = std::find_if(eclasses.cbegin(), eclasses.cend(), [](auto&& a1) { return a1.second.size() > 1; });
  }
  // eliminate duplicate conditions
  std::sort(conds.begin(), conds.end());
  conds.erase(std::unique(conds.begin(), conds.end()), conds.end());

  for (const auto& cond : conds) {
    QueryVertexID u = sorted_v_[cond.first];
    QueryVertexID v = sorted_v_[cond.second];
    result.emplace_back(u, v);
  }
  return result;
}

std::map<QueryVertexID, std::set<QueryVertexID>> AutomorphismCheck::getAEquivalenceClasses(
    const std::vector<std::vector<QueryVertexID>>& A) {
  std::map<QueryVertexID, std::set<QueryVertexID>> eclasses;
  for (unsigned int id = 0; id < bliss_q_.get_nof_vertices(); ++id) {
    std::set<QueryVertexID> eclass;
    eclass.insert(id);
    for (auto& perm : A) {
      eclass.insert(perm[id]);
    }
    QueryVertexID rep = *eclass.begin();  // the smallest element
    eclasses[rep].insert(eclass.cbegin(), eclass.cend());
  }
  return eclasses;
}

}  // namespace circinus

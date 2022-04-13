#include "algorithms/automorphism_check.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <unordered_set>
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

PartialOrder AutomorphismCheck::getPartialOrder() {
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
      DLOG(INFO) << "equivalence classes:" << ss.str();
    }

    QueryVertexID n0 = *eclass.cbegin();
    for (auto&& f : automorphisms) {
      QueryVertexID min = *std::min_element(std::next(eclass.cbegin()), eclass.cend(),
                                            [&f](QueryVertexID n, QueryVertexID m) { return f[n] < f[m]; });
      conds.emplace_back(n0, min);
      if (verbosePlannerLog()) {
        DLOG(INFO) << "add partial order condition " << n0 << '<' << min;
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
      DLOG(INFO) << "automorphisms:" << ss.str();
    }

    eclasses = getAEquivalenceClasses(automorphisms);
    eclass_it = std::find_if(eclasses.cbegin(), eclasses.cend(), [](auto&& a1) { return a1.second.size() > 1; });
  }
  // eliminate duplicate conditions
  std::sort(conds.begin(), conds.end());
  conds.erase(std::unique(conds.begin(), conds.end()), conds.end());

  return orderVertices(conds);
}

PartialOrder AutomorphismCheck::orderVertices(
    const std::vector<std::pair<QueryVertexID, QueryVertexID>>& sorted_conds) {
  auto n_qvs = q_->getNumVertices();
  PartialOrder po(n_qvs);

  if (sorted_conds.empty()) return po;

  auto& smaller_larger_qvs = po.constraints;  // qv: {smaller qvs, larger qvs}
  std::queue<QueryVertexID> queue;
  std::vector<uint32_t> checked_smaller_count(n_qvs, 0);  // i: # smaller qvs already checked regarding query vertices i
  unordered_set<QueryVertexID> covered_smaller_qvs;       // tmp
  uint32_t ordered_count = 0;

  /* find qvs with constraints */
  for (const auto& cond : sorted_conds) {
    QueryVertexID u = sorted_v_[cond.first];
    QueryVertexID v = sorted_v_[cond.second];
    smaller_larger_qvs[u].second.insert(v);
    smaller_larger_qvs[v].first.insert(u);
  }

  /* order qvs from smaller to larger by partial order */
  for (auto& p : smaller_larger_qvs) {  // find qvs without constraint of smaller qvs
    if (p.second.first.empty()) {
      // a query vertex is put in queue when all smaller qvs are poped and ordered
      queue.push(p.first);
    }
  }
  while (!queue.empty()) {
    auto current_qv = queue.front();
    queue.pop();
    auto& entry = smaller_larger_qvs.at(current_qv);
    covered_smaller_qvs.clear();
    // find all smaller qvs than `smaller` and deduplicate constraints
    for (auto smaller : entry.first) {
      // here `smaller_than_smaller` is an exhaustive set of smaller qvs regarding `smaller`
      auto& smaller_than_smaller = smaller_larger_qvs.at(smaller).first;
      covered_smaller_qvs.insert(smaller_than_smaller.begin(), smaller_than_smaller.end());
    }
    for (auto smaller : entry.first) {
      if (covered_smaller_qvs.insert(smaller).second) {  // not covered by any other smaller qv
        po.po_constraint_adj[current_qv].push_back(smaller);
      }
    }
    // update smaller qvs to an exhaustive set
    entry.first.swap(covered_smaller_qvs);
    po.order_index[current_qv] = ordered_count++;
    // check all larger qvs and update queue
    for (auto next_qv : entry.second) {
      if (++checked_smaller_count[next_qv] == smaller_larger_qvs.at(next_qv).first.size()) {
        queue.push(next_qv);
      }
    }
  }
  DCHECK_EQ(ordered_count, smaller_larger_qvs.size());

  // make exhaustive set of larger qvs
  for (auto& p : smaller_larger_qvs) {
    for (auto smaller : p.second.first) {
      smaller_larger_qvs.at(smaller).second.insert(p.first);
    }
  }
  for (auto& p : smaller_larger_qvs) {
    po.n_all_related_constraints[p.first] += p.second.first.size() + p.second.second.size();
  }
  return po;
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

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

#include "algorithms/partial_order.h"

#include <algorithm>
#include <iostream>
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"

namespace circinus {

std::ostream& PartialOrder::printMinimum(std::ostream& stream) const {
  for (QueryVertexID u = 0; u < po_constraint_adj.size(); ++u) {
    if (po_constraint_adj[u].empty()) continue;
    for (auto smaller : po_constraint_adj[u]) {
      stream << ' ' << smaller << '<' << u;
    }
  }
  return stream;
}

std::ostream& PartialOrder::printFull(std::ostream& stream) const {
  for (auto& p : constraints) {
    stream << "\nvertex " << p.first << ": " << n_all_related_constraints[p.first] << "\t<(";
    for (auto larger : p.second.second) {
      stream << larger << ' ';
    }
    stream << ")\t>(";
    for (auto smaller : p.second.first) {
      stream << smaller << ' ';
    }
    stream << ")";
  }
  return stream;
}

PartialOrder::EnforcePlan PartialOrder::orderVertices(const std::vector<QueryVertexID>& vertices) const {
  EnforcePlan plan;
  if (po_constraint_adj.empty()) return plan;
  auto n = vertices.size();
  auto tmp_order = vertices;
  std::vector<std::vector<uint32_t>> tmp_constraints_adj(n);

  std::sort(tmp_order.begin(), tmp_order.end(),
            [this](QueryVertexID a, QueryVertexID b) { return order_index[a] < order_index[b]; });

  unordered_map<QueryVertexID, uint32_t> seen;
  // some vertices may be unconstrained within the given set
  std::vector<bool> is_constrained(n, false);
  for (uint32_t i = 0; i < n; ++i) {
    auto current_qv = tmp_order[i];
    if (!seen.empty()) {
      for (auto v : po_constraint_adj[current_qv]) {
        auto pos = seen.find(v);
        if (pos != seen.end()) {
          is_constrained[pos->second] = true;
          tmp_constraints_adj[i].push_back(pos->second);
        }
      }
    }
    is_constrained[i] = !tmp_constraints_adj[i].empty();
    seen.emplace(current_qv, i);
  }

  std::vector<uint32_t> new_order(n);
  uint32_t n_constrained = 0;
  uint32_t n_unconstrained = 0;
  for (uint32_t i = 0; i < n; ++i) {
    if (is_constrained[i]) {
      new_order[i] = n_constrained;
      ++n_constrained;
    } else {
      new_order[i] = n_unconstrained + n;
      ++n_unconstrained;
    }
  }
  if (n_constrained == 0) {
    return plan;
  } else {
    plan.has_constraint = true;
  }
  plan.order.resize(n);
  plan.constraints_adj.resize(n);
  for (uint32_t i = 0; i < n; ++i) {
    if (new_order[i] > n) {
      new_order[i] -= n - n_constrained;
    }
    plan.order[new_order[i]] = tmp_order[i];
  }

  for (uint32_t i = 0; i < n; ++i) {
    if (!tmp_constraints_adj[i].empty()) {
      auto& adj = plan.constraints_adj[new_order[i]];
      adj.reserve(tmp_constraints_adj[i].size());
      for (auto index : tmp_constraints_adj[i]) {
        plan.constraints_adj[new_order[i]].push_back(new_order[index]);
      }
    }
  }
  DCHECK(checkEnforcePlan(plan)) << "bug in orderVertices";
  return plan;
}

bool PartialOrder::checkEnforcePlan(const EnforcePlan& plan) const {
  bool has_constraint = false;
  for (uint32_t i = 0; i < plan.order.size(); ++i) {
    auto v = plan.order[i];
    if (plan.constraints_adj[i].empty()) continue;
    has_constraint = true;
    auto& smaller_indices = plan.constraints_adj[i];
    auto pos = constraints.find(v);
    CHECK(pos != constraints.end());
    // not checking for missing constraint now
    for (auto index : smaller_indices) {
      CHECK_LT(index, i);                                       // indices of smaller vertices should be smaller
      CHECK_EQ(pos->second.first.count(plan.order[index]), 1);  // no non-existing constraint
    }
  }
  return plan.has_constraint == has_constraint;
}

std::vector<std::pair<bool, uint32_t>> PartialOrder::getConstraintsForVertex(
    QueryVertexID v, const unordered_map<QueryVertexID, uint32_t>& seen_vertices) const {
  return getConstraintsForVertexImpl<false>(v, seen_vertices);
}

std::vector<std::pair<bool, QueryVertexID>> PartialOrder::getConstraintVerticesForVertex(
    QueryVertexID v, const unordered_map<QueryVertexID, uint32_t>& seen_vertices) const {
  return getConstraintsForVertexImpl<true>(v, seen_vertices);
}

template <bool return_qv>
std::vector<std::pair<bool, std::conditional_t<return_qv, QueryVertexID, uint32_t>>>
PartialOrder::getConstraintsForVertexImpl(QueryVertexID v,
                                          const unordered_map<QueryVertexID, uint32_t>& seen_vertices) const {
  auto pos = constraints.find(v);
  if (pos == constraints.end()) return {};

  std::vector<std::pair<bool, std::conditional_t<return_qv, QueryVertexID, uint32_t>>> res;
  auto & [ smaller_vs, larger_vs ] = pos->second;
  // handle vertices smaller than v
  std::vector<QueryVertexID> relevant_vs;
  for (auto smaller : smaller_vs) {
    if (seen_vertices.count(smaller)) {
      relevant_vs.push_back(smaller);
    }
  }
  if (!relevant_vs.empty()) {
    std::sort(relevant_vs.begin(), relevant_vs.end(),
              [this](QueryVertexID a, QueryVertexID b) { return order_index[a] < order_index[b]; });
    for (uint32_t i = 0; i < relevant_vs.size(); ++i) {
      bool uncovered = true;
      for (uint32_t j = i + 1; j < relevant_vs.size(); ++j) {
        // check whether the i-th is already covered by (smaller than) the j-th vertex
        if (constraints.at(relevant_vs[j]).first.count(relevant_vs[i]) == 1) {
          uncovered = false;
          break;
        }
      }
      if (uncovered) {
        if
          constexpr(return_qv) {
            res.emplace_back(false, relevant_vs[i]);
            continue;
          }
        res.emplace_back(false, seen_vertices.at(relevant_vs[i]));
      }
    }
  }
  // handle vertices larger than v
  relevant_vs.clear();
  for (auto larger : larger_vs) {
    if (seen_vertices.count(larger)) {
      relevant_vs.push_back(larger);
    }
  }
  if (!relevant_vs.empty()) {
    // greatest comes first
    std::sort(relevant_vs.begin(), relevant_vs.end(),
              [this](QueryVertexID a, QueryVertexID b) { return order_index[a] > order_index[b]; });
    for (uint32_t i = 0; i < relevant_vs.size(); ++i) {
      bool uncovered = true;
      for (uint32_t j = i + 1; j < relevant_vs.size(); ++j) {
        // check whether the i-th is already covered by (larger than) the j-th vertex
        if (constraints.at(relevant_vs[j]).second.count(relevant_vs[i]) == 1) {
          uncovered = false;
          break;
        }
      }
      if (uncovered) {
        if
          constexpr(return_qv) {
            res.emplace_back(true, relevant_vs[i]);
            continue;
          }
        res.emplace_back(true, seen_vertices.at(relevant_vs[i]));
      }
    }
  }
  return res;
}

}  // namespace circinus

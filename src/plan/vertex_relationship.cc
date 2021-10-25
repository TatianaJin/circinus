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

#include "plan/vertex_relationship.h"

#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "utils/flags.h"
#include "utils/utils.h"

namespace circinus {

std::pair<QueryVertexID, std::vector<QueryVertexID>> VertexRelationship::findReusableSet(
    QueryVertexID target, std::vector<QueryVertexID>& set_vertices,
    const unordered_set<QueryVertexID>& existing_vertices) {
  DLOG(INFO) << "findReusableSet for " << target << " among [" << toString(set_vertices) << " ]";
  std::pair<QueryVertexID, std::vector<QueryVertexID>> res;
  res.first = DUMMY_QUERY_VERTEX;
  std::vector<QueryVertexID>& uncovered_parent_vertices = res.second;
  auto q = qv_equivalence_->getQueryGraph();
  auto po = qv_equivalence_->getPartialOrder();

  auto target_nbrs = q->getOutNeighbors(target);
  for (uint32_t i = 0; i < target_nbrs.second; ++i) {
    if (existing_vertices.count(target_nbrs.first[i])) {
      uncovered_parent_vertices.push_back(target_nbrs.first[i]);
    }
  }
  std::vector<QueryVertexID> enforced_target_constraints;
  auto& target_constraints = po->po_constraint_adj[target];
  for (uint32_t i = 0; i < target_constraints.size(); ++i) {
    if (existing_vertices.count(target_constraints[i])) {
      enforced_target_constraints.push_back(target_constraints[i]);
    }
  }

  // For now consider only equivalence but not containment
  for (QueryVertexID set_vertex : set_vertices) {
    if (q->getVertexLabel(set_vertex) != q->getVertexLabel(target)) continue;
    bool is_equivalent = true;
    // neighborhood equivalence
    auto set_nbrs = q->getOutNeighbors(set_vertex);
    uint32_t target_parent_index = 0;
    for (uint32_t i = 0; i < set_nbrs.second; ++i) {
      if (existing_vertices.count(set_nbrs.first[i])) {
        if (target_parent_index < uncovered_parent_vertices.size() &&
            uncovered_parent_vertices[target_parent_index] == set_nbrs.first[i]) {
          ++target_parent_index;
        } else {
          is_equivalent = false;
          break;
        }
      }
    }
    // partial order constraint equivalence
    if (is_equivalent) {
      auto& set_constraints = po->po_constraint_adj[set_vertex];
      uint32_t target_constraint_index = 0;
      for (uint32_t i = 0; i < set_constraints.size(); ++i) {
        if (existing_vertices.count(set_constraints[i])) {
          if (target_constraint_index < enforced_target_constraints.size() &&
              enforced_target_constraints[target_constraint_index] == set_constraints[i]) {
            ++target_constraint_index;
          } else {
            is_equivalent = false;
            break;
          }
        }
      }
    }
    // if not check target_parent_index == uncovered_parent_vertices.size(),
    // target neighbors can be a superset of reusable set neighbors
    if (is_equivalent && target_parent_index == uncovered_parent_vertices.size()) {
      res.first = set_vertex;
      uncovered_parent_vertices.clear();
      if (verbosePlannerLog()) {
        LOG(INFO) << "target " << target << " reusable set vertex " << set_vertex;
      }
      return res;
    }
  }
  uncovered_parent_vertices.clear();
  return res;
}

}  // namespace circinus

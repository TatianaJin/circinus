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
    const unordered_set<QueryVertexID>& existing_vertices) const {
  DLOG(INFO) << "findReusableSet for " << target << " among [" << toString(set_vertices) << " ]";
  std::pair<QueryVertexID, std::vector<QueryVertexID>> res;
  res.first = DUMMY_QUERY_VERTEX;
  auto q = qv_equivalence_->getQueryGraph();
  auto po = qv_equivalence_->getPartialOrder();

  unordered_set<QueryVertexID> uncovered_parent_vertices;
  auto target_nbrs = q->getOutNeighbors(target);
  for (uint32_t i = 0; i < target_nbrs.second; ++i) {
    if (existing_vertices.count(target_nbrs.first[i])) {
      uncovered_parent_vertices.insert(target_nbrs.first[i]);
    }
  }
  auto target_constraints = po->constraints.find(target);

  for (QueryVertexID set_vertex : set_vertices) {
    if (q->getVertexLabel(set_vertex) != q->getVertexLabel(target)) continue;
    bool is_reusable = true;
    LOG(INFO) << "check neighbor of " << target << " and " << set_vertex;
    // neighborhood equivalence
    auto set_nbrs = q->getOutNeighbors(set_vertex);
    if (set_nbrs.second > target_nbrs.second) continue;  // the degree filter of set_vertex will make it unusable
    auto uncovered = uncovered_parent_vertices;
    for (uint32_t i = 0; i < set_nbrs.second; ++i) {
      if (existing_vertices.count(set_nbrs.first[i])) {
        auto pos = uncovered.find(set_nbrs.first[i]);
        if (pos != uncovered.end()) {
          uncovered.erase(pos);
        } else {
          is_reusable = false;
          break;
        }
      }
    }
    if (uncovered.size() == uncovered_parent_vertices.size()) continue;
    if (!is_reusable) continue;

    LOG(INFO) << "check po gt of " << target << " and " << set_vertex;
    // partial order constraint equivalence
    auto set_constraints = po->constraints.find(set_vertex);
    if (set_constraints != po->constraints.end()) {
      for (auto cond : set_constraints->second.first) {
        if (existing_vertices.count(cond) == 0) continue;
        if (target_constraints == po->constraints.end() || target_constraints->second.first.count(cond) == 0) {
          is_reusable = false;
          break;
        }
      }
      if (!is_reusable) continue;
      LOG(INFO) << "check po lt of " << target << " and " << set_vertex;
      for (auto cond : set_constraints->second.second) {
        if (existing_vertices.count(cond) == 0) continue;
        if (target_constraints == po->constraints.end() || target_constraints->second.second.count(cond) == 0) {
          is_reusable = false;
          break;
        }
      }
    }

    // if not check target_parent_index == uncovered_parent_vertices.size(),
    // target neighbors can be a superset of reusable set neighbors
    if (is_reusable) {
      res.first = set_vertex;
      res.second.insert(res.second.end(), uncovered.begin(), uncovered.end());
      if (verbosePlannerLog()) {
        LOG(INFO) << "target " << target << " reusable set vertex " << set_vertex << " uncovered parents "
                  << toString(res.second);
      }
      return res;
    }
  }
  uncovered_parent_vertices.clear();
  return res;
}

}  // namespace circinus

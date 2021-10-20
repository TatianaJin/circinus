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

  auto target_nbrs = q->getOutNeighbors(target);
  for (uint32_t i = 0; i < target_nbrs.second; ++i) {
    if (existing_vertices.count(target_nbrs.first[i])) {
      uncovered_parent_vertices.push_back(target_nbrs.first[i]);
    }
  }

  // For now consider only equivalence but not containment
  for (QueryVertexID set_vertex : set_vertices) {
    if (q->getVertexLabel(set_vertex) != q->getVertexLabel(target)) continue;
    auto set_nbrs = q->getOutNeighbors(set_vertex);
    uint32_t target_parent_index = 0;
    bool is_equivalent = true;
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

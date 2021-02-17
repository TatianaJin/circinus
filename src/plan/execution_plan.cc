// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "plan/execution_plan.h"

namespace circinus {

void ExecutionPlan::populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                                         const std::vector<int>& cover_table) {
  query_graph_ = g;
  matching_order_ = matching_order;
  cover_table_ = cover_table;
  operators_.reserve(g->getNumVertices());

  auto parent = matching_order.front();
  if (matching_order.size() == 1) {
    // TODO(tatiana): handle no traversal
    return;
  }
  // for traversal, we expand to one vertex at a time according to matching_order
  auto root = newExpandEdgeOperator(parent, matching_order[1]);
  Operator *current, *prev = root;
  std::unordered_set<QueryVertexID> existing_vertices;
  existing_vertices.insert(parent);
  existing_vertices.insert(matching_order[1]);

  std::vector<QueryVertexID> parents;
  parents.reserve(g->getNumVertices() - 1);
  for (uint32_t i = 2; i < matching_order.size(); ++i) {
    auto target_vertex = matching_order[i];
    // find parent vertices
    auto neighbors = g->getOutNeighbors(target_vertex);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (existing_vertices.count(neighbors.first[j]) != 0) {
        parents.push_back(neighbors.first[j]);
      }
    }
    if (parents.size() == 1) {  // only one parent, ExpandEdge
      current = newExpandEdgeOperator(parents.front(), target_vertex);
    } else {  // more than one parents, ExpandVertex (use set intersection)
      current = newExpandVertexOperator(parents, target_vertex);
    }
    prev->setNext(current);
    prev = current;
    parents.clear();
    existing_vertices.insert(target_vertex);
  }
  DCHECK_EQ(existing_vertices.size(), g->getNumVertices());
}

ExpandEdgeOperator* ExecutionPlan::newExpandEdgeOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex) {
  auto ret = new ExpandEdgeOperator(query_graph_, parent_vertex, target_vertex, cover_table_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

ExpandVertexOperator* ExecutionPlan::newExpandVertexOperator(std::vector<QueryVertexID>& parents,
                                                             QueryVertexID target_vertex) {
  auto ret = new ExpandVertexOperator(query_graph_, parents, target_vertex, cover_table_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

}  // namespace circinus

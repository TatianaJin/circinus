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

#include "ops/operators.h"

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

  /* for traversal, we expand to one vertex at a time according to matching_order */
  uint32_t n_keys = 0, n_sets = 0;
  Operator *current, *prev = nullptr;
  // existing vertices in current query subgraph
  std::unordered_set<QueryVertexID> existing_vertices;
  std::vector<QueryVertexID> parents;
  parents.reserve(g->getNumVertices() - 1);

  // handle first vertex
  query_vertex_indices_[parent] = 0;
  n_keys += (cover_table[parent] == 1);
  n_sets += (cover_table[parent] != 1);
  existing_vertices.insert(parent);
  // handle following traversals
  for (uint32_t i = 1; i < matching_order.size(); ++i) {
    auto target_vertex = matching_order[i];
    // record query vertex index
    if (cover_table[target_vertex] == 1) {
      query_vertex_indices_[target_vertex] = n_keys;
      ++n_keys;
    } else {
      query_vertex_indices_[target_vertex] = n_sets;
      ++n_sets;
    }

    // find parent vertices
    if (i == 1) {  // the first edge
      prev = newExpandEdgeOperator(parent, target_vertex);
    } else {
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
    }

    existing_vertices.insert(target_vertex);
  }

  DCHECK_EQ(existing_vertices.size(), g->getNumVertices());
}

TraverseOperator* ExecutionPlan::newExpandEdgeOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex) {
  auto ret =
      ExpandEdgeOperator::newExpandEdgeOperator(parent_vertex, target_vertex, cover_table_, query_vertex_indices_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandVertexOperator(std::vector<QueryVertexID>& parents,
                                                         QueryVertexID target_vertex) {
  auto ret = ExpandVertexOperator::newExpandVertexOperator(query_graph_, parents, target_vertex, cover_table_,
                                                           query_vertex_indices_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

}  // namespace circinus

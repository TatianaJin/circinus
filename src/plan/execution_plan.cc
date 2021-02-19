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
#include "ops/expand_into_operator.h"
#include "ops/operators.h"
#include "ops/expand_key_key_vertex_operator.h"
#include "ops/expand_set_vertex_operator.h"
#include "ops/expand_set_to_key_vertex_operator.h"

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
  uint32_t n_keys = 0, n_sets = 0;
  Operator *current, *prev = nullptr;
  std::unordered_set<QueryVertexID> existing_vertices;
  existing_vertices.insert(parent);
  existing_vertices.insert(matching_order[1]);

  std::vector<QueryVertexID> key_parents;
  std::vector<QueryVertexID> set_parents;
  key_parents.reserve(g->getNumVertices() - 1);
  set_parents.reserve(g->getNumVertices() - 1);
  for (uint32_t i = 2; i < matching_order.size(); ++i) {
    auto target_vertex = matching_order[i];
    // record query vertex index
    if (cover_table[target_vertex] == 1) {
      query_vertex_indices_[target_vertex] = n_keys;
      ++n_keys;
    } else {
      query_vertex_indices_[target_vertex] = n_sets;
      ++n_sets;
    }

    auto neighbors = g->getOutNeighbors(target_vertex);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (existing_vertices.count(neighbors.first[j]) != 0) {
        if (cover_table[neighbors.first[j]]) {
          key_parents.push_back(neighbors.first[j]);
        } else {
          set_parents.push_back(neighbors.first[j]);
        }
      }
    }
    
    if (key_parents.size() + set_parents.size() == 1) {  // only one parent, ExpandEdge
      auto front = key_parents.size() == 1 ? key_parents.front() : set_parents.front();
      current = newExpandEdgeOperator(front, target_vertex);
      prev->setNext(current);
      prev = current;
    } else {  // more than one parents, ExpandVertex (use set intersection)
      if (cover_table[target_vertex]) {
        if (key_parents.size() != 0 && set_parents.size() != 0) { 
          current = newExpandKeyKeyVertexOperator(key_parents, target_vertex);
          prev->setNext(current);
          prev = current;
          current = newExpandIntoOperator(set_parents, target_vertex);
          prev->setNext(current);
          prev = current;
        } else if (key_parents.size() != 0) {
          current = newExpandKeyKeyVertexOperator(key_parents, target_vertex);
          prev->setNext(current);
          prev = current;
        } else {
          current = newExpandSetToKeyVertexOperator(key_parents, target_vertex);
          prev->setNext(current);
          prev = current;
        }
      } else {
        current = newExpandSetVertexOperator(key_parents, target_vertex);
        prev->setNext(current);
        prev = current;
      }
    }
    key_parents.clear();
    set_parents.clear();
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

TraverseOperator* ExecutionPlan::newExpandIntoOperator(std::vector<QueryVertexID>& parents,
                                                             QueryVertexID target_vertex) {
  auto ret = new ExpandIntoOperator(query_graph_, parents, target_vertex, cover_table_, query_vertex_indices_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandSetVertexOperator(std::vector<QueryVertexID>& parents,
                                                             QueryVertexID target_vertex) {
  auto ret = new ExpandSetVertexOperator(query_graph_, parents, target_vertex, cover_table_, query_vertex_indices_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandSetToKeyVertexOperator(std::vector<QueryVertexID>& parents,
                                                             QueryVertexID target_vertex) {
  auto ret = new ExpandSetToKeyVertexOperator(query_graph_, parents, target_vertex, cover_table_, query_vertex_indices_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandKeyKeyVertexOperator(std::vector<QueryVertexID>& parents,
                                                             QueryVertexID target_vertex) {
  auto ret = new ExpandKeyKeyVertexOperator(query_graph_, parents, target_vertex, cover_table_, query_vertex_indices_);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

}  // namespace circinus

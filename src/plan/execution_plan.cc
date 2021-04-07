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

#include <array>
#include <unordered_set>
#include <vector>

#include "ops/operators.h"
#include "utils/hashmap.h"

namespace circinus {

void ExecutionPlan::populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                                         const std::vector<int>& cover_table, Profiler* profiler) {
  query_graph_ = g;
  cover_table_ = cover_table;
  operators_.reserve(g->getNumVertices());
  operators_.setProfiler(profiler);
  root_query_vertex_ = matching_order.front();
  if (matching_order.size() == 1) {
    // TODO(tatiana): handle no traversal
    return;
  }

  /* for traversal, we expand to one vertex at a time according to matching_order */
  uint32_t n_keys = 0, n_sets = 0;
  Operator *current, *prev = nullptr;
  // existing vertices in current query subgraph
  unordered_set<QueryVertexID> existing_vertices;
  // label: {set index}, {key index}
  unordered_map<LabelID, std::array<std::vector<uint32_t>, 2>> label_existing_vertices_indices;
  std::array<std::vector<QueryVertexID>, 2> parents;
  auto& key_parents = parents[1];
  auto& set_parents = parents[0];
  key_parents.reserve(g->getNumVertices() - 1);
  set_parents.reserve(g->getNumVertices() - 1);

  // handle first vertex
  auto parent = matching_order.front();
  query_vertex_indices_[parent] = 0;
  n_keys += (cover_table[parent] == 1);
  n_sets += (cover_table[parent] != 1);
  existing_vertices.insert(parent);
  label_existing_vertices_indices[query_graph_->getVertexLabel(parent)][cover_table[parent] == 1].push_back(0);

  // handle following traversals
  for (uint32_t i = 1; i < matching_order.size(); ++i) {
    auto target_vertex = matching_order[i];
    auto target_label = query_graph_->getVertexLabel(target_vertex);
    const auto& same_label_v_indices = label_existing_vertices_indices[target_label];
    // record query vertex index
    if (cover_table[target_vertex] == 1) {
      query_vertex_indices_[target_vertex] = n_keys;
      ++n_keys;
    } else {
      query_vertex_indices_[target_vertex] = n_sets;
      ++n_sets;
    }

    if (i == 1) {  // the first edge: target has one and only one parent
      prev = newExpandEdgeOperator(parent, target_vertex, cover_table, same_label_v_indices);
    } else {
      // find parent vertices
      auto neighbors = g->getOutNeighbors(target_vertex);
      for (uint32_t j = 0; j < neighbors.second; ++j) {
        if (existing_vertices.count(neighbors.first[j]) != 0) {
          parents[cover_table[neighbors.first[j]] == 1].push_back(neighbors.first[j]);
        }
      }
      // create operators
      if (key_parents.size() + set_parents.size() == 1) {  // only one parent, ExpandEdge
        auto front = key_parents.size() == 1 ? key_parents.front() : set_parents.front();
        current = newExpandEdgeOperator(front, target_vertex, cover_table, same_label_v_indices);
      } else {  // more than one parents, ExpandVertex (use set intersection)
        if (cover_table[target_vertex] == 1) {
          if (key_parents.size() != 0 && set_parents.size() != 0) {
            if (key_parents.size() == 1) {
              current = newExpandEdgeOperator(key_parents.front(), target_vertex, cover_table, same_label_v_indices);
            } else {
              current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_v_indices);
            }
            prev->setNext(current);
            prev = current;
            current = newExpandIntoOperator(set_parents, target_vertex, key_parents);
          } else if (key_parents.size() != 0) {
            current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_v_indices);
          } else {
            current = newExpandSetToKeyVertexOperator(set_parents, target_vertex, same_label_v_indices);
          }
        } else {
          current = newExpandSetVertexOperator(key_parents, target_vertex, same_label_v_indices);
        }
      }
      prev->setNext(current);
      prev = current;
      key_parents.clear();
      set_parents.clear();
    }

    existing_vertices.insert(target_vertex);
    label_existing_vertices_indices[target_label][cover_table[target_vertex] == 1].push_back(
        query_vertex_indices_[target_vertex]);
  }

  // output
  auto output_op = newOutputOperator();
  prev->setNext(output_op);

  DCHECK_EQ(existing_vertices.size(), g->getNumVertices());
}

void ExecutionPlan::populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                                         const std::vector<int>& cover_table,
                                         const unordered_map<QueryVertexID, uint32_t>& level_become_key) {
  query_graph_ = g;
  cover_table_ = cover_table;
  operators_.reserve(g->getNumVertices());
  root_query_vertex_ = matching_order.front();
  dynamic_cover_key_level_ = level_become_key;
  if (matching_order.size() == 1) {
    // TODO(tatiana): handle no traversal
    return;
  }

  uint32_t n_keys = 0, n_sets = 0;
  std::vector<QueryVertexID> set_vertices;
  set_vertices.reserve(matching_order.size() - dynamic_cover_key_level_.size());
  Operator *current, *prev = nullptr;
  // existing vertices in current query subgraph
  unordered_set<QueryVertexID> existing_vertices;
  unordered_map<LabelID, std::vector<uint32_t>> label_existing_vertices_map;
  std::array<std::vector<QueryVertexID>, 2> parents;
  auto& key_parents = parents[1];
  auto& set_parents = parents[0];
  key_parents.reserve(g->getNumVertices() - 1);
  set_parents.reserve(g->getNumVertices() - 1);

  // handle first vertex
  auto parent = matching_order.front();
  query_vertex_indices_[parent] = 0;
  if (dynamic_cover_key_level_.count(parent) == 0) {
    set_vertices.push_back(parent);
    ++n_sets;
  } else {
    ++n_keys;
  }
  existing_vertices.insert(parent);
  label_existing_vertices_map[query_graph_->getVertexLabel(parent)].push_back(parent);
  std::vector<std::vector<QueryVertexID>> add_keys_at_level(matching_order.size());
  unordered_map<QueryVertexID, uint32_t> input_query_vertex_indices;

  // handle following traversals
  for (uint32_t i = 1; i < matching_order.size(); ++i) {
    // FIXME(tatiana): for debug: check set vertices and query_vertex_indices_ are in correspondence
    for (uint32_t set_i = 0; set_i < set_vertices.size(); ++set_i) {
      CHECK_EQ(query_vertex_indices_[set_vertices[set_i]], set_i);
    }
    auto target_vertex = matching_order[i];
    auto target_label = query_graph_->getVertexLabel(target_vertex);

    const auto& same_label_vertices = label_existing_vertices_map[target_label];
    std::array<std::vector<uint32_t>, 2> same_label_indices;
    for (auto v : same_label_vertices) {
      // note that here cover_table_ reflects the compression key of the input instead of the output
      same_label_indices[cover_table_[v] == 1].push_back(query_vertex_indices_[v]);
    }

    // record query vertex index
    auto key_level_pos = dynamic_cover_key_level_.find(target_vertex);
    if (key_level_pos == dynamic_cover_key_level_.end()) {  // target vertex is in set
      DCHECK_NE(cover_table_[target_vertex], 1);
      if (!add_keys_at_level[i].empty()) {
        input_query_vertex_indices = query_vertex_indices_;
        for (auto v : add_keys_at_level[i]) {
          auto& v_idx = query_vertex_indices_[v];
          CHECK_EQ(v, set_vertices[v_idx]);
          // swap with the last set vertex
          set_vertices[v_idx] = set_vertices.back();
          set_vertices.pop_back();
          query_vertex_indices_[set_vertices[v_idx]] = v_idx;
          v_idx = n_keys++;
          cover_table_[v] = 1;
        }
        n_sets -= add_keys_at_level[i].size();
      }
      query_vertex_indices_[target_vertex] = n_sets;
      ++n_sets;
      set_vertices.push_back(target_vertex);
    } else if (key_level_pos->second <= i) {  // target vertex is in key
      query_vertex_indices_[target_vertex] = n_keys;
      ++n_keys;
    } else {  // target vertex is in key in some later subqueries
      /* now we assume the existing vertices will not be turned into key if target vertex is in the vertex cover */
      CHECK(add_keys_at_level[i].empty());
      query_vertex_indices_[target_vertex] = n_sets;
      ++n_sets;
      set_vertices.push_back(target_vertex);
      add_keys_at_level[key_level_pos->second].push_back(target_vertex);  // add to key later
      cover_table_[target_vertex] = 0;
    }

    if (i == 1) {  // the first edge: target has one and only one parent
      prev = newExpandEdgeOperator(parent, target_vertex, cover_table_, same_label_indices);
    } else {
      // find parent vertices
      auto neighbors = g->getOutNeighbors(target_vertex);
      for (uint32_t j = 0; j < neighbors.second; ++j) {
        if (existing_vertices.count(neighbors.first[j])) {
          parents[cover_table_[neighbors.first[j]] == 1].push_back(neighbors.first[j]);
        }
      }
      // create operators
      if (cover_table_[target_vertex] == 1) {  // target is in key
        CHECK(add_keys_at_level[i].empty());
        if (key_parents.size() != 0 && set_parents.size() != 0) {
          if (key_parents.size() == 1) {
            current = newExpandEdgeOperator(key_parents.front(), target_vertex, cover_table_, same_label_indices);
          } else {
            current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_indices);
          }
          prev->setNext(current);
          prev = current;
          current = newExpandIntoOperator(set_parents, target_vertex, key_parents);
        } else if (key_parents.size() == 1) {
          current = newExpandEdgeOperator(key_parents.front(), target_vertex, cover_table_, same_label_indices);
        } else if (key_parents.size() > 1) {
          current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_indices);
        } else if (set_parents.size() == 1) {
          current = newExpandEdgeOperator(set_parents.front(), target_vertex, cover_table_, same_label_indices);
        } else {
          current = newExpandSetToKeyVertexOperator(set_parents, target_vertex, same_label_indices);
        }
      } else {  // target is in set, then  all parents should be in key, and key enumeration may be needed
        if (!add_keys_at_level[i].empty()) {  // key enumeration is needed
          current = newEnumerateKeyExpandToSetOperator(key_parents, target_vertex, add_keys_at_level[i],
                                                       input_query_vertex_indices, same_label_indices);
          add_keys_at_level[i].clear();
        } else if (key_parents.size() == 1) {
          current = newExpandEdgeOperator(key_parents.front(), target_vertex, cover_table_, same_label_indices);
        } else {
          current = newExpandSetVertexOperator(key_parents, target_vertex, same_label_indices);
        }
      }
      prev->setNext(current);
      prev = current;
      key_parents.clear();
      set_parents.clear();
    }
    existing_vertices.insert(target_vertex);
    label_existing_vertices_map[target_label].push_back(target_vertex);
  }

  // output
  auto output_op = newOutputOperator();
  prev->setNext(output_op);

  DCHECK_EQ(existing_vertices.size(), g->getNumVertices());
}

TraverseOperator* ExecutionPlan::newExpandEdgeOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                       const std::vector<int>& cover_table,
                                                       const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  auto ret = ExpandEdgeOperator::newExpandEdgeOperator(parent_vertex, target_vertex, cover_table, query_vertex_indices_,
                                                       same_label_indices[1], same_label_indices[0],
                                                       getSetPruningThreshold(target_vertex));
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandIntoOperator(const std::vector<QueryVertexID>& parents,
                                                       QueryVertexID target_vertex,
                                                       const std::vector<QueryVertexID>& prev_key_parents) {
  auto ret = new ExpandIntoOperator(parents, target_vertex, query_vertex_indices_, prev_key_parents);
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandSetVertexOperator(
    std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  auto ret = new ExpandKeyToSetVertexOperator(parents, target_vertex, query_vertex_indices_, same_label_indices[1],
                                              same_label_indices[0], getSetPruningThreshold(target_vertex));
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandSetToKeyVertexOperator(
    std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  auto ret = new ExpandSetToKeyVertexOperator(parents, target_vertex, query_vertex_indices_, same_label_indices[1],
                                              same_label_indices[0], getSetPruningThreshold(target_vertex));
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandKeyKeyVertexOperator(
    std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  auto ret = new ExpandKeyToKeyVertexOperator(parents, target_vertex, query_vertex_indices_, same_label_indices[1],
                                              same_label_indices[0], getSetPruningThreshold(target_vertex));
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

Operator* ExecutionPlan::newOutputOperator() {
  auto ret = OutputOperator::newOutputOperator(OutputType::Count, &outputs_);
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newEnumerateKeyExpandToSetOperator(
    const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::vector<QueryVertexID>& keys_to_enumerate,
    unordered_map<QueryVertexID, uint32_t> input_query_vertex_indices,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  auto ret = new EnumerateKeyExpandToSetOperator(
      parents, target_vertex, input_query_vertex_indices, query_vertex_indices_, keys_to_enumerate, cover_table_,
      same_label_indices[1], same_label_indices[0], getSetPruningThreshold(target_vertex));
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

}  // namespace circinus

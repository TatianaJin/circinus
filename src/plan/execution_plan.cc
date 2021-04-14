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

#include <algorithm>
#include <array>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/bipartite_graph.h"
#include "ops/operators.h"
#include "utils/hashmap.h"

namespace circinus {

void addBipartiteGraphToOperator(
    QueryVertexID qv1, QueryVertexID qv2, TraverseOperator* op) {
  op->addBipartiteGraph(new BipartiteGraph(qv1, qv2));
}

// BipartiteGraph used in this function but not the other function with the same name
void ExecutionPlan::populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                                         const std::vector<int>& cover_table, Profiler* profiler) {
  matching_order_ = matching_order;
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
      if (cover_table[target_vertex] != 1) {
        prev = newExpandEdgeKeyToSetOperator(parent, target_vertex, same_label_v_indices);
      } else if (cover_table[parent] == 1) {
        prev = newExpandEdgeKeyToKeyOperator(parent, target_vertex, same_label_v_indices);
      } else {
        prev = newExpandEdgeSetToKeyOperator(parent, target_vertex, same_label_v_indices, std::vector<uint32_t>{});
      }
      addBipartiteGraphToOperator(parent, target_vertex, dynamic_cast<TraverseOperator*>(prev));
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
        if (cover_table[target_vertex] != 1) {
          current = newExpandEdgeKeyToSetOperator(key_parents.front(), target_vertex, same_label_v_indices);
        } else if (key_parents.size() == 1) {
          current = newExpandEdgeKeyToKeyOperator(key_parents.front(), target_vertex, same_label_v_indices);
        } else {
          current = newExpandEdgeSetToKeyOperator(
              set_parents.front(), target_vertex, same_label_v_indices,
              label_existing_vertices_indices[query_graph_->getVertexLabel(set_parents.front())][0]);
        }
      } else {  // more than one parents, ExpandVertex (use set intersection)
        if (cover_table[target_vertex] == 1) {
          if (key_parents.size() != 0 && set_parents.size() != 0) {
            if (key_parents.size() == 1) {
              current = newExpandEdgeKeyToKeyOperator(key_parents.front(), target_vertex, same_label_v_indices);
            } else {
              current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_v_indices);
            }
            prev->setNext(current);
            prev = current;
            current = newExpandIntoOperator(set_parents, target_vertex, key_parents, label_existing_vertices_indices);
          } else if (key_parents.size() != 0) {
            current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_v_indices);
          } else {
            current = newExpandSetToKeyVertexOperator(set_parents, target_vertex, same_label_v_indices,
                                                      label_existing_vertices_indices);
          }
        } else {
          current = newExpandSetVertexOperator(key_parents, target_vertex, same_label_v_indices);
        }
      }
      // careful with this key-first-set-later order
      for (auto& qv : key_parents) {
        addBipartiteGraphToOperator(qv, target_vertex, dynamic_cast<TraverseOperator*>(current));
      }
      for (auto& qv : set_parents) {
        addBipartiteGraphToOperator(qv, target_vertex, dynamic_cast<TraverseOperator*>(current));
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
  std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> same_label_indices;  // {{keys},{sets}}
  for (auto& pair : label_existing_vertices_indices) {
    if ((!pair.second[1].empty() && !pair.second[0].empty()) || pair.second[0].size() > 1) {
      same_label_indices.emplace_back(std::move(pair.second[1]), std::move(pair.second[0]));
    }
  }
  auto output_op = newOutputOperator(std::move(same_label_indices));
  prev->setNext(output_op);

  DCHECK_EQ(existing_vertices.size(), g->getNumVertices());
}

void ExecutionPlan::populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                                         const std::vector<int>& cover_table,
                                         const unordered_map<QueryVertexID, uint32_t>& level_become_key) {
  matching_order_ = matching_order;
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
  set_vertices.reserve(matching_order.size());
  Operator *current, *prev = nullptr;
  // existing vertices in current query subgraph
  unordered_set<QueryVertexID> existing_vertices;
  unordered_map<LabelID, std::vector<QueryVertexID>> label_existing_vertices_map;
  std::array<std::vector<QueryVertexID>, 2> parents;
  auto& key_parents = parents[1];
  auto& set_parents = parents[0];
  key_parents.reserve(g->getNumVertices() - 1);
  set_parents.reserve(g->getNumVertices() - 1);
  std::vector<std::vector<QueryVertexID>> add_keys_at_level(matching_order.size());

  // handle first vertex
  auto parent = matching_order.front();
  query_vertex_indices_[parent] = 0;
  if (dynamic_cover_key_level_.count(parent) == 0) {
    set_vertices.push_back(parent);
    ++n_sets;
  } else {
    if (dynamic_cover_key_level_.find(parent)->second == 0) {
      ++n_keys;
    } else {
      ++n_sets;
      set_vertices.push_back(parent);
      add_keys_at_level[dynamic_cover_key_level_.find(parent)->second].push_back(parent);  // add to key later
      cover_table_[parent] = 0;
    }
  }

  existing_vertices.insert(parent);
  label_existing_vertices_map[query_graph_->getVertexLabel(parent)].push_back(parent);
  unordered_map<QueryVertexID, uint32_t> input_query_vertex_indices;

  // handle following traversals
  for (uint32_t i = 1; i < matching_order.size(); ++i) {
    for (uint32_t set_i = 0; set_i < set_vertices.size(); ++set_i) {
      DCHECK_EQ(query_vertex_indices_[set_vertices[set_i]], set_i);
    }
    auto target_vertex = matching_order[i];
    auto target_label = query_graph_->getVertexLabel(target_vertex);

    const auto& same_label_vertices = label_existing_vertices_map[target_label];
    std::array<std::vector<uint32_t>, 2> same_label_indices;  // indices in the input
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
        addKeys(add_keys_at_level[i], set_vertices, n_keys);
        n_sets -= add_keys_at_level[i].size();
      }
      query_vertex_indices_[target_vertex] = n_sets;
      ++n_sets;
      set_vertices.push_back(target_vertex);
    } else if (key_level_pos->second <= i) {  // target vertex is in key
      query_vertex_indices_[target_vertex] = n_keys;
      ++n_keys;
    } else {  // target vertex is in key in some later subqueries
      if (!add_keys_at_level[i].empty()) {
        input_query_vertex_indices = query_vertex_indices_;
        addKeys(add_keys_at_level[i], set_vertices, n_keys);
        n_sets -= add_keys_at_level[i].size();
      }
      query_vertex_indices_[target_vertex] = n_sets;
      ++n_sets;
      set_vertices.push_back(target_vertex);
      add_keys_at_level[key_level_pos->second].push_back(target_vertex);  // add to key later
      cover_table_[target_vertex] = 0;
    }

    // find parent vertices
    auto neighbors = g->getOutNeighbors(target_vertex);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (existing_vertices.count(neighbors.first[j])) {
        parents[cover_table_[neighbors.first[j]] == 1].push_back(neighbors.first[j]);
      }
    }
    // print cover for current subquery
    std::stringstream cover_ss;
    for (uint32_t k = 0; k <= i; ++k) {
      cover_ss << ' ' << matching_order_[k] << ':' << cover_table_[matching_order_[k]];
    }
    LOG(INFO) << (i + 1) << "-cover>" << cover_ss.str();

    if (i == 1) {  // the first edge: target has one and only one parent
      if (cover_table_[target_vertex] != 1 && !add_keys_at_level[i].empty()) {
        prev = newEnumerateKeyExpandToSetOperator(key_parents, target_vertex, add_keys_at_level[i],
                                                  input_query_vertex_indices, same_label_indices,
                                                  label_existing_vertices_map);
        add_keys_at_level[i].clear();
      } else if (cover_table_[target_vertex] != 1) {
        prev = newExpandEdgeKeyToSetOperator(parent, target_vertex, same_label_indices);
      } else if (cover_table_[parent] == 1) {
        prev = newExpandEdgeKeyToKeyOperator(parent, target_vertex, same_label_indices);
      } else {
        prev = newExpandEdgeSetToKeyOperator(parent, target_vertex, same_label_indices, std::vector<uint32_t>{});
      }
    } else {
      // create operators
      if (cover_table_[target_vertex] == 1) {  // target is in key
        CHECK(add_keys_at_level[i].empty()) << i << " size " << add_keys_at_level[i].size();
        if (key_parents.size() != 0 && set_parents.size() != 0) {
          if (key_parents.size() == 1) {
            current = newExpandEdgeKeyToKeyOperator(key_parents.front(), target_vertex, same_label_indices);
          } else {
            current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_indices);
          }
          prev->setNext(current);
          prev = current;
          current = newExpandIntoOperator(set_parents, target_vertex, key_parents, label_existing_vertices_map);
        } else if (key_parents.size() == 1) {
          current = newExpandEdgeKeyToKeyOperator(key_parents.front(), target_vertex, same_label_indices);
        } else if (key_parents.size() > 1) {
          current = newExpandKeyKeyVertexOperator(key_parents, target_vertex, same_label_indices);
        } else if (set_parents.size() == 1) {
          current = newExpandEdgeSetToKeyOperator(set_parents.front(), target_vertex, same_label_indices,
                                                  label_existing_vertices_map);
        } else {
          current = newExpandSetToKeyVertexOperator(set_parents, target_vertex, same_label_indices,
                                                    label_existing_vertices_map);
        }
      } else {  // target is in set, then  all parents should be in key, and key enumeration may be needed
        if (!add_keys_at_level[i].empty()) {  // key enumeration is needed
          current = newEnumerateKeyExpandToSetOperator(key_parents, target_vertex, add_keys_at_level[i],
                                                       input_query_vertex_indices, same_label_indices,
                                                       label_existing_vertices_map);
          add_keys_at_level[i].clear();
        } else if (key_parents.size() == 1) {
          current = newExpandEdgeKeyToSetOperator(key_parents.front(), target_vertex, same_label_indices);
        } else {
          current = newExpandSetVertexOperator(key_parents, target_vertex, same_label_indices);
        }
      }
      prev->setNext(current);
      prev = current;
    }
    key_parents.clear();
    set_parents.clear();
    existing_vertices.insert(target_vertex);
    label_existing_vertices_map[target_label].push_back(target_vertex);
  }

  // output
  std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> same_label_indices;  // {{keys},{sets}}
  for (auto& pair : label_existing_vertices_map) {
    std::vector<uint32_t> same_label_sets;
    std::vector<uint32_t> same_label_keys;
    for (auto v : pair.second) {
      if (cover_table_[v] == 1) {
        same_label_keys.push_back(query_vertex_indices_[v]);
      } else {  // set vertex
        same_label_sets.push_back(query_vertex_indices_[v]);
      }
    }
    if ((!same_label_keys.empty() && !same_label_sets.empty()) || same_label_sets.size() > 1) {
      same_label_indices.emplace_back(std::move(same_label_keys), std::move(same_label_sets));
    }
  }
  auto output_op = newOutputOperator(std::move(same_label_indices));
  prev->setNext(output_op);

  DCHECK_EQ(existing_vertices.size(), g->getNumVertices());
}

TraverseOperator* ExecutionPlan::newExpandEdgeKeyToKeyOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  SubgraphFilter* filter = nullptr;
#ifdef USE_FILTER
  filter = SubgraphFilter::newSetPrunningSubgraphFilter(same_label_indices[0]);
  subgraph_filters_.push_back(filter);
#endif
  auto ret = ExpandEdgeOperator::newExpandEdgeKeyToKeyOperator(parent_vertex, target_vertex, query_vertex_indices_,
                                                               same_label_indices[1], same_label_indices[0],
                                                               getSetPruningThreshold(target_vertex), filter);
  setMatchingOrderIndices(target_vertex, ret);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandEdgeKeyToSetOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  SubgraphFilter* filter = nullptr;
#ifdef USE_FILTER
  if (same_label_indices[0].empty()) {
    subgraph_filters_.push_back(SubgraphFilter::newDummyFilter());
  } else {
    std::vector<std::vector<uint32_t>> pruning_sets(1);
    pruning_sets.front().resize(same_label_indices[0].size() + 1);
    std::copy(same_label_indices[0].begin(), same_label_indices[0].end(), pruning_sets.front().begin());
    pruning_sets.front().back() = query_vertex_indices_[target_vertex];
    subgraph_filters_.push_back(SubgraphFilter::newSetPrunningSubgraphFilter(std::move(pruning_sets)));
  }
  filter = subgraph_filters_.back();
#endif
  auto ret = ExpandEdgeOperator::newExpandEdgeKeyToSetOperator(parent_vertex, target_vertex, query_vertex_indices_,
                                                               same_label_indices[1], same_label_indices[0],
                                                               getSetPruningThreshold(target_vertex), filter);
  setMatchingOrderIndices(target_vertex, ret);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandEdgeSetToKeyOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& target_same_label_indices,
    const std::vector<uint32_t>& parent_same_label_set_indices) {
  std::vector<uint32_t> target_same_label_set_indices = target_same_label_indices[0];
  // if parent is in the same label indices, remove
  if (query_graph_->getVertexLabel(parent_vertex) == query_graph_->getVertexLabel(target_vertex)) {
    auto parent_index = query_vertex_indices_[parent_vertex];
    for (uint32_t i = 0; i < target_same_label_set_indices.size(); ++i) {
      if (target_same_label_set_indices[i] == parent_index) {
        target_same_label_set_indices[i] = target_same_label_set_indices.back();
        break;
      }
    }
    target_same_label_set_indices.pop_back();
  }
#ifdef USE_FILTER
  if (target_same_label_indices[0].size() < 2 && parent_same_label_set_indices.size() < 2) {
    subgraph_filters_.push_back(SubgraphFilter::newDummyFilter());
  } else if (target_same_label_indices[0].size() < 2) {
    subgraph_filters_.push_back(SubgraphFilter::newSetPrunningSubgraphFilter(parent_same_label_set_indices));
  } else if (parent_same_label_set_indices.size() < 2) {
    subgraph_filters_.push_back(SubgraphFilter::newSetPrunningSubgraphFilter(target_same_label_indices[0]));
  } else {
    // if parent and target are of the same label, only one group of pruning sets
    if (query_graph_->getVertexLabel(parent_vertex) == query_graph_->getVertexLabel(target_vertex)) {
      subgraph_filters_.push_back(SubgraphFilter::newSetPrunningSubgraphFilter(target_same_label_indices[0]));
    } else {
      std::vector<std::vector<uint32_t>> pruning_sets{target_same_label_indices[0], parent_same_label_set_indices};
      subgraph_filters_.push_back(SubgraphFilter::newSetPrunningSubgraphFilter(std::move(pruning_sets)));
    }
  }
  auto ret = ExpandEdgeOperator::newExpandEdgeSetToKeyOperator(
      parent_vertex, target_vertex, query_vertex_indices_, target_same_label_indices[1], target_same_label_set_indices,
      getSetPruningThreshold(target_vertex), subgraph_filters_.back());
#else
  auto ret = ExpandEdgeOperator::newExpandEdgeSetToKeyOperator(
      parent_vertex, target_vertex, query_vertex_indices_, target_same_label_indices[1], target_same_label_set_indices,
      getSetPruningThreshold(target_vertex), nullptr);
#endif
  setMatchingOrderIndices(target_vertex, ret);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandIntoOperator(const std::vector<QueryVertexID>& parents,
                                                       QueryVertexID target_vertex,
                                                       const std::vector<QueryVertexID>& prev_key_parents,
                                                       std::vector<std::vector<uint32_t>>&& pruning_set_indices) {
  SubgraphFilter* filter = nullptr;
#ifdef USE_FILTER
  filter = createFilter(std::move(pruning_set_indices));
#endif
  auto ret = new ExpandIntoOperator(parents, target_vertex, query_vertex_indices_, prev_key_parents, filter);
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandSetVertexOperator(
    std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  SubgraphFilter* filter = nullptr;
#ifdef USE_FILTER
  if (same_label_indices[0].empty()) {
    filter = SubgraphFilter::newDummyFilter();
  } else {
    std::vector<std::vector<uint32_t>> pruning_sets(1);
    pruning_sets.front().resize(same_label_indices[0].size() + 1);
    std::copy(same_label_indices[0].begin(), same_label_indices[0].end(), pruning_sets.front().begin());
    pruning_sets.front().back() = query_vertex_indices_[target_vertex];
    filter = SubgraphFilter::newSetPrunningSubgraphFilter(std::move(pruning_sets));
  }
  subgraph_filters_.push_back(filter);
#endif
  auto ret = new ExpandKeyToSetVertexOperator(parents, target_vertex, query_vertex_indices_, same_label_indices[1],
                                              same_label_indices[0], getSetPruningThreshold(target_vertex), filter);
  setMatchingOrderIndices(target_vertex, ret);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandSetToKeyVertexOperator(
    const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices,
    std::vector<std::vector<uint32_t>>&& pruning_set_indices) {
  std::vector<uint32_t> target_same_label_set_indices = same_label_indices[0];

  auto target_label = query_graph_->getVertexLabel(target_vertex);
  unordered_set<uint32_t> same_label_parent_indices;
  for (auto parent : parents) {
    if (query_graph_->getVertexLabel(parent) == target_label) {
      same_label_parent_indices.insert(query_vertex_indices_[parent]);
    }
  }
  // if parent is in the same label indices, remove
  if (!same_label_parent_indices.empty()) {
    target_same_label_set_indices.erase(
        std::remove_if(
            target_same_label_set_indices.begin(), target_same_label_set_indices.end(),
            [&same_label_parent_indices](uint32_t idx) { return same_label_parent_indices.count(idx) == 1; }),
        target_same_label_set_indices.end());
  }
  SubgraphFilter* filter = nullptr;
#ifdef USE_FILTER
  filter = createFilter(std::move(pruning_set_indices));
#endif
  auto ret =
      new ExpandSetToKeyVertexOperator(parents, target_vertex, query_vertex_indices_, same_label_indices[1],
                                       target_same_label_set_indices, getSetPruningThreshold(target_vertex), filter);
  setMatchingOrderIndices(target_vertex, ret);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newExpandKeyKeyVertexOperator(
    std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices) {
  SubgraphFilter* filter = nullptr;
#ifdef USE_FILTER
  filter = SubgraphFilter::newSetPrunningSubgraphFilter(same_label_indices[0]);
  subgraph_filters_.push_back(filter);
#endif
  auto ret = new ExpandKeyToKeyVertexOperator(parents, target_vertex, query_vertex_indices_, same_label_indices[1],
                                              same_label_indices[0], getSetPruningThreshold(target_vertex), filter);
  setMatchingOrderIndices(target_vertex, ret);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

Operator* ExecutionPlan::newOutputOperator(
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>&& same_label_indices) {
  auto ret = OutputOperator::newOutputOperator(OutputType::Count, &outputs_, std::move(same_label_indices));
  operators_.push_back(ret);
  return ret;
}

TraverseOperator* ExecutionPlan::newEnumerateKeyExpandToSetOperator(
    const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const std::vector<QueryVertexID>& keys_to_enumerate,
    unordered_map<QueryVertexID, uint32_t> input_query_vertex_indices,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices,
    const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map) {
  unordered_map<LabelID, int> pruning_labels;              // label: index of pruning sets
  std::vector<std::vector<uint32_t>> pruning_set_indices;  // indices in the output
  // add indices of sets whose labels are the same with target_vertex
  LabelID target_label = query_graph_->getVertexLabel(target_vertex);
  auto insert = pruning_labels.insert({target_label, -1});
  if (same_label_indices[0].size() - keys_to_enumerate.size() > 0) {
    std::vector<uint32_t> target_same_label_set_indices;  // indices in the output
    DCHECK_EQ(label_existing_vertices_map.count(target_label), 1);
    for (auto v : label_existing_vertices_map.at(target_label)) {
      if (cover_table_[v] != 1) {
        target_same_label_set_indices.push_back(query_vertex_indices_.at(v));
      }
    }
#ifdef USE_FILTER
    target_same_label_set_indices.push_back(query_vertex_indices_.at(target_vertex));
#endif
    pruning_set_indices.emplace_back(std::move(target_same_label_set_indices));
    insert.first->second = 0;
  }
  auto enumerated_key_pruning_indices =
      getPruningSets(keys_to_enumerate, label_existing_vertices_map, pruning_labels, pruning_set_indices);
  std::vector<uint64_t> pruning_set_thresholds(pruning_set_indices.size(), FLAGS_set_pruning_threshold);
  SubgraphFilter* filter = createFilter(std::move(pruning_set_indices));
  if (FLAGS_set_pruning_threshold == 0) {
    for (auto& pair : pruning_labels) {
      if (pair.second != -1) {
        pruning_set_thresholds[pair.second] = query_graph_->getVertexCardinalityByLabel(pair.first);
      }
    }
  }
  filter->setPruningSetThresholds(std::move(pruning_set_thresholds));
  auto ret = new EnumerateKeyExpandToSetOperator(parents, target_vertex, input_query_vertex_indices,
                                                 query_vertex_indices_, keys_to_enumerate, cover_table_,
                                                 same_label_indices, std::move(enumerated_key_pruning_indices), filter);
  setMatchingOrderIndices(target_vertex, ret);
  target_vertex_to_ops_[target_vertex] = ret;
  operators_.push_back(ret);
  return ret;
}

}  // namespace circinus

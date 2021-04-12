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

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/operator.h"
#include "ops/output_operator.h"
#include "ops/traverse_operator.h"
#include "plan/operator_tree.h"
#include "utils/flags.h"
#include "utils/hashmap.h"
#include "utils/profiler.h"

namespace circinus {

class ExecutionPlan {
 private:
  const QueryGraph* query_graph_;

  std::vector<std::vector<VertexID>> candidate_sets_;
  unordered_map<QueryVertexID, Operator*> target_vertex_to_ops_;
  OperatorTree operators_;
  Outputs outputs_;
  std::vector<SubgraphFilter*> subgraph_filters_;  // owned, need to delete upon destruction

  QueryVertexID root_query_vertex_;
  std::vector<int> cover_table_;
  unordered_map<QueryVertexID, uint32_t> dynamic_cover_key_level_;

  /** The index of each query vertex in the CompressedSubgraphs
   * For key vertices, the index is n_keys-th key following the matching order
   * For non-key vertices, the index n_sets-th key following the matching order */
  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;

  void addKeys(const std::vector<QueryVertexID>& keys_to_add, std::vector<QueryVertexID>& set_vertices,
               uint32_t& n_keys) {
    for (auto v : keys_to_add) {
      auto& v_idx = query_vertex_indices_[v];
      // swap with the last set vertex
      set_vertices[v_idx] = set_vertices.back();
      set_vertices.pop_back();
      query_vertex_indices_[set_vertices[v_idx]] = v_idx;
      v_idx = n_keys++;
      cover_table_[v] = 1;
    }
  }

 public:
  ~ExecutionPlan() {
    for (auto filter : subgraph_filters_) {
      delete filter;
    }
    subgraph_filters_.clear();
  }

  void populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                            const std::vector<int>& cover_table, Profiler* profiler = nullptr);

  void populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                            const std::vector<int>& cover_table,
                            const unordered_map<QueryVertexID, uint32_t>& level_become_key);

  void setProfiler(Profiler* profiler) { operators_.setProfiler(profiler); }

  void printPhysicalPlan() const {
    if (operators_.empty()) return;
    printOperatorChain(operators_.root());
  }

  void printProfiledPlan(std::ostream& oss) const {
    if (operators_.empty()) return;
    uint32_t op_idx = 0;
    auto op = operators_.root();
    while (op != nullptr) {
      oss << op_idx << "," << op->toProfileString() << std::endl;
      ++op_idx;
      op = op->getNext();
    }
  }

  void printLabelFrequency() const {
    unordered_map<LabelID, std::vector<QueryVertexID>> labels;
    uint32_t n_keys = 0;
    for (QueryVertexID v = 0; v < query_graph_->getNumVertices(); ++v) {
      labels[query_graph_->getVertexLabel(v)].push_back(v);
      n_keys += isInCover(v);
    }
    for (auto& pair : labels) {
      std::stringstream ss;
      for (auto v : pair.second) {
        ss << ' ' << v << '(' << (isInCover(v) ? "key" : "set") << ')';
      }
      LOG(INFO) << "label " << pair.first << ':' << ss.str();
    }
    LOG(INFO) << n_keys << " keys " << (query_graph_->getNumVertices() - n_keys) << " sets";
  }

  inline const OperatorTree& getOperators() const { return operators_; }
  inline OperatorTree& getOperators() { return operators_; }

  inline void setCandidateSets(std::vector<std::vector<VertexID>>& cs) {
    candidate_sets_ = cs;
    uint64_t total_key_candidate_size = 1, total_set_candidate_size = 1;
    for (uint32_t i = 0; i < candidate_sets_.size(); ++i) {
      DLOG(INFO) << "query vertex " << i << ": " << candidate_sets_[i].size() << " candidates";
      if (isInCover(i)) {
        total_key_candidate_size *= candidate_sets_[i].size();
      } else {
        total_set_candidate_size *= candidate_sets_[i].size();
      }
      if (i == root_query_vertex_) continue;
      ((TraverseOperator*)target_vertex_to_ops_[i])->setCandidateSets(&candidate_sets_[i]);
    }
    LOG(INFO) << "total_key_candidate_size=" << total_key_candidate_size
              << ",total_set_candidate_size=" << total_set_candidate_size;
    LOG(INFO) << "ratio " << ((double)total_set_candidate_size / total_key_candidate_size);
  }

  inline const std::vector<VertexID>& getCandidateSet(QueryVertexID id) const {
    DCHECK_LT(id, candidate_sets_.size());
    return candidate_sets_[id];
  }
  
  inline const std::vector<std::vector<VertexID>>& getCandidateSets() const {
    return candidate_sets_;
  }


  inline const bool isToKey(QueryVertexID id) const { return dynamic_cover_key_level_.count(id) != 0; }

  inline const uint32_t getToKeyLevel(QueryVertexID id) const { return dynamic_cover_key_level_.find(id)->second; }

  inline OperatorTree cloneOperators() const { return operators_.clone(); }

  inline QueryVertexID getRootQueryVertexID() const { return root_query_vertex_; }

  inline bool isInCover(QueryVertexID id) const { return cover_table_[id] == 1; }

  inline const Outputs& getOutputs() const { return outputs_; }
  inline Outputs& getOutputs() { return outputs_; }

 protected:
  inline SubgraphFilter* createFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices) {
    SubgraphFilter* filter = nullptr;
    if (pruning_set_indices.empty()) {
      filter = SubgraphFilter::newDummyFilter();
    } else {
      filter = SubgraphFilter::newSetPrunningSubgraphFilter(std::move(pruning_set_indices));
    }
    subgraph_filters_.push_back(filter);
    return filter;
  }

  inline void getPruningSetsForParents(
      const std::vector<QueryVertexID>& parents,
      const unordered_map<LabelID, std::array<std::vector<uint32_t>, 2>>& label_existing_vertices_indices,
      unordered_set<LabelID>& pruning_labels, std::vector<std::vector<uint32_t>>& pruning_set_indices) {
    for (auto& p : parents) {
      auto label = query_graph_->getVertexLabel(p);
      if (pruning_labels.insert(label).second) {
        auto pos = label_existing_vertices_indices.find(label);
        if (pos != label_existing_vertices_indices.end() && pos->second[0].size() > 1) {
          pruning_set_indices.push_back(pos->second[0]);
        }
      }
    }
  }

  inline void getPruningSetsForParents(const std::vector<QueryVertexID>& parents,
                                       const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map,
                                       unordered_set<LabelID>& pruning_labels,
                                       std::vector<std::vector<uint32_t>>& pruning_set_indices) const {
    for (auto& p : parents) {
      auto label = query_graph_->getVertexLabel(p);
      if (pruning_labels.insert(label).second) {
        auto pos = label_existing_vertices_map.find(label);
        if (pos != label_existing_vertices_map.end() && pos->second.size() > 1) {
          std::vector<uint32_t> indices;
          for (auto v : pos->second) {
            // note that here cover_table_ reflect the cover in the target-expanded output
            if (cover_table_[v] != 1) {
              indices.push_back(query_vertex_indices_.at(v));
            }
          }
          if (indices.size() > 1) {
            pruning_set_indices.push_back(std::move(indices));
          }
        }
      }
    }
  }

  inline std::vector<int> getPruningSets(
      const std::vector<QueryVertexID>& query_vertices,
      const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map,
      unordered_map<LabelID, int>& pruning_labels, std::vector<std::vector<uint32_t>>& pruning_set_indices) const {
    std::vector<int> vertex_pruning_set_indices;
    vertex_pruning_set_indices.resize(query_vertices.size(), -1);
    for (uint32_t i = 0; i < query_vertices.size(); ++i) {
      auto label = query_graph_->getVertexLabel(query_vertices[i]);
      auto inserted = pruning_labels.insert({label, -1});
      if (inserted.second) {
        auto pos = label_existing_vertices_map.find(label);
        if (pos != label_existing_vertices_map.end() && pos->second.size() > 1) {
          std::vector<uint32_t> indices;
          for (auto v : pos->second) {
            // note that here cover_table_ reflect the cover in the target-expanded output
            if (cover_table_[v] != 1) {
              indices.push_back(query_vertex_indices_.at(v));
            }
          }
          if (indices.size() > 1) {
            vertex_pruning_set_indices[i] = pruning_set_indices.size();
            inserted.first->second = pruning_set_indices.size();
            pruning_set_indices.push_back(std::move(indices));
          }
        }
      } else {
        vertex_pruning_set_indices[i] = inserted.first->second;
      }
    }
    return vertex_pruning_set_indices;
  }

  inline uint64_t getSetPruningThreshold(QueryVertexID pruning_qv) {
    return FLAGS_set_pruning_threshold == 0
               ? query_graph_->getVertexCardinalityByLabel(query_graph_->getVertexLabel(pruning_qv))
               : FLAGS_set_pruning_threshold;
  }

  inline void printOperatorChain(const Operator* root, const std::string& indent = "") const {
    LOG(INFO) << indent << root->toString();
    if (root->getNext() != nullptr) {
      printOperatorChain(root->getNext(), indent + "| ");
    }
  }

  TraverseOperator* newExpandEdgeKeyToKeyOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                  const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandEdgeKeyToSetOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                  const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandEdgeSetToKeyOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                  const std::array<std::vector<uint32_t>, 2>& target_same_label_indices,
                                                  const std::vector<uint32_t>& parent_same_label_indices);
  TraverseOperator* newExpandKeyKeyVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                                  const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandSetVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                               const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandSetToKeyVertexOperator(const std::vector<QueryVertexID>& parents,
                                                    QueryVertexID target_vertex,
                                                    const std::array<std::vector<uint32_t>, 2>& same_label_indices,
                                                    std::vector<std::vector<uint32_t>>&& pruning_set_indices);
  TraverseOperator* newExpandIntoOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                          const std::vector<QueryVertexID>& prev_key_parents,
                                          std::vector<std::vector<uint32_t>>&& pruning_set_indices);
  TraverseOperator* newEnumerateKeyExpandToSetOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::vector<QueryVertexID>& keys_to_enumerate,
      unordered_map<QueryVertexID, uint32_t> input_query_vertex_indices,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices,
      const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map);
  Operator* newOutputOperator(std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>&&);

  inline TraverseOperator* newExpandEdgeSetToKeyOperator(
      QueryVertexID parent_vertex, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& target_same_label_indices,
      const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map) {
    std::vector<uint32_t> parent_prune;
#ifdef USE_FILTER
    auto parent_label = query_graph_->getVertexLabel(parent_vertex);
    auto parent_pos = label_existing_vertices_map.find(parent_label);
    if (parent_label != query_graph_->getVertexLabel(target_vertex) &&
        parent_pos != label_existing_vertices_map.end() && parent_pos->second.size() > 1) {
      for (auto v : parent_pos->second) {
        if (cover_table_[v] != 1) {
          parent_prune.push_back(query_vertex_indices_[v]);
        }
      }
    }
#endif
    return newExpandEdgeSetToKeyOperator(parent_vertex, target_vertex, target_same_label_indices, parent_prune);
  }

  inline TraverseOperator* newExpandSetToKeyVertexOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices,
      const unordered_map<LabelID, std::array<std::vector<uint32_t>, 2>>& label_existing_vertices_indices) {
#ifdef USE_FILTER
    unordered_set<LabelID> pruning_labels;
    std::vector<std::vector<uint32_t>> pruning_set_indices;
    pruning_labels.insert(query_graph_->getVertexLabel(target_vertex));
    if (same_label_indices[0].size() > 1) {
      pruning_set_indices.push_back(same_label_indices[0]);
    }
    getPruningSetsForParents(parents, label_existing_vertices_indices, pruning_labels, pruning_set_indices);
    return newExpandSetToKeyVertexOperator(parents, target_vertex, same_label_indices, std::move(pruning_set_indices));
#else
    return newExpandSetToKeyVertexOperator(parents, target_vertex, same_label_indices,
                                           std::vector<std::vector<uint32_t>>{});
#endif
  }

  inline TraverseOperator* newExpandSetToKeyVertexOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices,
      const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map) {
#ifdef USE_FILTER
    unordered_set<LabelID> pruning_labels;
    std::vector<std::vector<uint32_t>> pruning_set_indices;
    pruning_labels.insert(query_graph_->getVertexLabel(target_vertex));
    if (same_label_indices[0].size() > 1) {
      pruning_set_indices.push_back(same_label_indices[0]);
    }
    getPruningSetsForParents(parents, label_existing_vertices_map, pruning_labels, pruning_set_indices);
    return newExpandSetToKeyVertexOperator(parents, target_vertex, same_label_indices, std::move(pruning_set_indices));
#else
    return newExpandSetToKeyVertexOperator(parents, target_vertex, same_label_indices,
                                           std::vector<std::vector<uint32_t>>{});
#endif
  }

  inline TraverseOperator* newExpandIntoOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::vector<QueryVertexID>& prev_key_parents,
      const unordered_map<LabelID, std::array<std::vector<uint32_t>, 2>>& label_existing_vertices_indices) {
#ifdef USE_FILTER
    unordered_set<LabelID> pruning_labels;
    std::vector<std::vector<uint32_t>> pruning_set_indices;
    getPruningSetsForParents(parents, label_existing_vertices_indices, pruning_labels, pruning_set_indices);
    return newExpandIntoOperator(parents, target_vertex, prev_key_parents, std::move(pruning_set_indices));
#else
    return newExpandIntoOperator(parents, target_vertex, prev_key_parents, std::vector<std::vector<uint32_t>>{});
#endif
  }

  inline TraverseOperator* newExpandIntoOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::vector<QueryVertexID>& prev_key_parents,
      const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map) {
#ifdef USE_FILTER
    unordered_set<LabelID> pruning_labels;
    std::vector<std::vector<uint32_t>> pruning_set_indices;
    getPruningSetsForParents(parents, label_existing_vertices_map, pruning_labels, pruning_set_indices);
    return newExpandIntoOperator(parents, target_vertex, prev_key_parents, std::move(pruning_set_indices));
#else
    return newExpandIntoOperator(parents, target_vertex, prev_key_parents, std::vector<std::vector<uint32_t>>{});
#endif
  }
};

}  // namespace circinus

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

#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/vertex_equivalence.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/operator.h"
#include "ops/output_operator.h"
#include "ops/traverse_operator.h"
#include "plan/vertex_relationship.h"
#include "utils/flags.h"
#include "utils/hashmap.h"

namespace circinus {

/** Uses flag:
 *    label_filter  Whether to set target and parent query vertex label for pruning adjacency lists before intersection.
 */
class ExecutionPlan {
 protected:
  const GraphType graph_type_;
  const QueryGraph* query_graph_;

  bool seperate_enumerate_ = FLAGS_seperate_enumeration;
  std::vector<Operator*> operators_;  // TODO(tatiana): use unique_ptr?

  /* matching order related */
  std::vector<QueryVertexID> matching_order_;

  /* compression related */
  bool inputs_are_keys_ = true;
  std::vector<int> cover_table_;
  unordered_map<QueryVertexID, uint32_t> dynamic_cover_key_level_;
  /** The index of each query vertex in the CompressedSubgraphs
   * For key vertices, the index is n_keys-th key following the matching order
   * For non-key vertices, the index n_sets-th key following the matching order */
  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;

  VertexRelationship* qv_relationship_ = nullptr;  // owned by NaivePlanner

  std::vector<double> step_costs_;

  /* partition query vertex's candidate search space */
  uint32_t partition_id_ = 0;  // TODO(tatiana): obsolete
  uint64_t query_cover_bits_ = 0;

  // transient variable for populating plan
  bool last_op_ = false;

  // for profile, subquery size, cover bits, remove parent edge flag
  std::vector<std::tuple<uint32_t, uint64_t, bool>> op_input_subquery_cover_;

  void transformSameLabelIndices(std::array<std::vector<uint32_t>, 2>& same_label_indices,
                                 unordered_set<QueryVertexID> keys_to_enumerate_set, uint32_t from) {
    auto it = same_label_indices[from].begin();
    while (it != same_label_indices[from].end()) {
      if (keys_to_enumerate_set.count(*it)) {
        same_label_indices[from ^ 1].push_back(*it);
        it = same_label_indices[from].erase(it);
      } else {
        ++it;
      }
    }
  }

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

  inline uint64_t getCoverBit(uint32_t size) const {
    uint64_t bits = 0;
    for (uint32_t i = 0; i < size; ++i) {
      if (isInCover(matching_order_[i])) {
        bits |= 1 << matching_order_[i];
      }
    }
    return bits;
  }

 public:
  explicit ExecutionPlan(GraphType graph_type = GraphType::Normal) : graph_type_(graph_type) {}

  ~ExecutionPlan() {
    for (auto& op : operators_) {
      delete op;
    }
    operators_.clear();
  }

  void populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                            const std::vector<int>& cover_table);

  void populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                            const std::vector<int>& cover_table,
                            const unordered_map<QueryVertexID, uint32_t>& level_become_key);

  inline void setQueryCoverBits(uint64_t bits) { query_cover_bits_ = bits; }
  inline uint64_t getQueryCoverBits() const { return query_cover_bits_; }
  void setInputAreKeys(bool flag) { inputs_are_keys_ = flag; }
  void setPartitionId(uint32_t partition_id) { partition_id_ = partition_id; }
  bool inputAreKeys() const { return inputs_are_keys_; }

  inline void setStepCosts(std::vector<double>&& step_costs) { step_costs_ = std::move(step_costs); }
  inline const auto& getStepCosts() const { return step_costs_; }

  inline void setVertexRelationship(VertexRelationship* vr) { qv_relationship_ = vr; }

  const auto& getOpInputSubqueryCover() const { return op_input_subquery_cover_; }

  void printPhysicalPlan() const {
    if (operators_.empty()) return;
    printOperatorChain(operators_[0]);
  }

  void printPhysicalPlan(std::ostream& oss) const {
    if (operators_.empty()) return;
    uint32_t op_idx = 0;
    auto op = operators_.front();
    while (op != nullptr) {
      oss << ++op_idx << "," << op->toString() << std::endl;
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

  inline uint32_t getPartitionId() const { return partition_id_; }

  inline const std::vector<Operator*>& getOperators() const { return operators_; }
  inline std::vector<Operator*>& getOperators() { return operators_; }

  inline QueryVertexID getRootQueryVertexID() const { return matching_order_.front(); }
  inline const auto& getMatchingOrder() const { return matching_order_; }

  inline bool isInCover(QueryVertexID id) const { return cover_table_[id] == 1; }
  inline const uint32_t getToKeyLevel(QueryVertexID id) const { return dynamic_cover_key_level_.find(id)->second; }

  inline uint32_t getQueryVertexOutputIndex(QueryVertexID qv) const { return query_vertex_indices_.at(qv); }

 protected:
  inline void setMatchingOrderIndices(QueryVertexID target_vertex, TraverseOperator* op) {
    std::vector<std::pair<bool, uint32_t>> matching_order_indices;
    for (auto& v : matching_order_) {
      matching_order_indices.emplace_back(cover_table_[v] == 1, query_vertex_indices_[v]);
      if (v == target_vertex) break;
    }
    op->setMatchingOrderIndices(std::move(matching_order_indices));
    op->setTargetLabel(FLAGS_label_filter ? query_graph_->getVertexLabel(target_vertex) : ALL_LABEL);
    op->setTargetDegree(query_graph_->getVertexOutDegree(target_vertex));
  }

  inline std::unique_ptr<SubgraphFilter> createFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices) {
    std::unique_ptr<SubgraphFilter> filter = nullptr;
    if (pruning_set_indices.empty()) {
      filter = SubgraphFilter::newDummyFilter();
    } else {
      filter = SubgraphFilter::newSetPrunningSubgraphFilter(std::move(pruning_set_indices));
    }
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
        CHECK(pos != label_existing_vertices_map.end());
        if (pos->second.size() > 1) {
          std::vector<uint32_t> indices;
          for (auto v : pos->second) {
            // note that here cover_table_ reflect the cover in the target-expanded output
            if (cover_table_[v] != 1) {
              indices.push_back(query_vertex_indices_.at(v));
            }
          }
          vertex_pruning_set_indices[i] = pruning_set_indices.size();
          inserted.first->second = pruning_set_indices.size();
          pruning_set_indices.push_back(std::move(indices));
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

  inline bool opToIntersectCandidates(QueryVertexID target_vertex) const {
    if (FLAGS_candidate_set_intersection == 2) {
      return dynamic_cover_key_level_.count(target_vertex) == 1;
    }
    if (FLAGS_candidate_set_intersection == 1) {
      return !last_op_;
    }
    if (FLAGS_candidate_set_intersection == 3) {
      return false;
    }
    CHECK_EQ(FLAGS_candidate_set_intersection, 0);
    return true;
  }

  virtual TraverseOperator* newExpandEdgeKeyToKeyOperator(
      QueryVertexID parent_vertex, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  virtual TraverseOperator* newExpandEdgeKeyToSetOperator(
      QueryVertexID parent_vertex, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices,
      unordered_map<QueryVertexID, uint32_t>& query_vertex_indices);
  virtual TraverseOperator* newExpandEdgeSetToKeyOperator(
      QueryVertexID parent_vertex, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& target_same_label_indices,
      const std::vector<uint32_t>& parent_same_label_indices);
  virtual TraverseOperator* newExpandKeyKeyVertexOperator(
      std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  virtual TraverseOperator* newExpandSetVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                                       const std::array<std::vector<uint32_t>, 2>& same_label_indices,
                                                       unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                                                       bool have_existing_target_set = false);
  virtual TraverseOperator* newExpandSetToKeyVertexOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices,
      std::vector<std::vector<uint32_t>>&& pruning_set_indices);
  virtual TraverseOperator* newExpandIntoOperator(const std::vector<QueryVertexID>& parents,
                                                  QueryVertexID target_vertex,
                                                  const std::vector<QueryVertexID>& prev_key_parents,
                                                  std::vector<std::vector<uint32_t>>&& pruning_set_indices);
  virtual TraverseOperator* newEnumerateKeyExpandToSetOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::vector<QueryVertexID>& keys_to_enumerate,
      unordered_map<QueryVertexID, uint32_t> input_query_vertex_indices,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices,
      const unordered_map<LabelID, std::vector<uint32_t>>& label_existing_vertices_map);

  virtual std::vector<TraverseOperator*> newExpandKeyToSetEnumerateKeyExpandToSetOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::vector<QueryVertexID>& keys_to_enumerate,
      unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
      std::array<std::vector<uint32_t>, 2>& same_label_indices,
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

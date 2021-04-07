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
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
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

  QueryVertexID root_query_vertex_;
  std::vector<int> cover_table_;
  unordered_map<QueryVertexID, uint32_t> dynamic_cover_key_level_;

  /** The index of each query vertex in the CompressedSubgraphs
   * For key vertices, the index is n_keys-th key following the matching order
   * For non-key vertices, the index n_sets-th key following the matching order */
  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;

 public:
  ~ExecutionPlan() {}

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
    candidate_sets_.swap(cs);
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

  inline OperatorTree cloneOperators() const { return operators_.clone(); }

  inline QueryVertexID getRootQueryVertexID() const { return root_query_vertex_; }

  inline bool isInCover(QueryVertexID id) const { return cover_table_[id] == 1; }

  inline const Outputs& getOutputs() const { return outputs_; }
  inline Outputs& getOutputs() { return outputs_; }

 private:
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

  TraverseOperator* newExpandEdgeOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                          const std::vector<int>& cover_table,
                                          const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandKeyKeyVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                                  const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandSetVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                               const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandSetToKeyVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                                    const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  TraverseOperator* newExpandIntoOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                          const std::vector<QueryVertexID>& prev_key_parents);
  TraverseOperator* newEnumerateKeyExpandToSetOperator(
      const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
      const std::vector<QueryVertexID>& keys_to_enumerate,
      unordered_map<QueryVertexID, uint32_t> input_query_vertex_indices,
      const std::array<std::vector<uint32_t>, 2>& same_label_indices);
  Operator* newOutputOperator();
};

}  // namespace circinus

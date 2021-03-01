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
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/operator.h"
#include "ops/traverse_operator.h"
#include "plan/operator_tree.h"

namespace circinus {

class ExecutionPlan {
 private:
  const QueryGraph* query_graph_;

  std::vector<std::vector<VertexID>> candidate_sets_;
  std::unordered_map<QueryVertexID, Operator*> target_vertex_to_ops_;
  OperatorTree operators_;

  QueryVertexID root_query_vertex_;
  std::vector<int> cover_table_;

  /** The index of each query vertex in the CompressedSubgraphs
   * For key vertices, the index is n_keys-th key following the matching order
   * For non-key vertices, the index n_sets-th key following the matching order */
  std::unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;

 public:
  ~ExecutionPlan() {}

  void populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                            const std::vector<int>& cover_table);

  void printPhysicalPlan() const {
    if (operators_.empty()) return;
    printOperatorChain(operators_.root());
  }

  inline const OperatorTree& getOperators() const { return operators_; }

  inline void setCandidateSets(std::vector<std::vector<VertexID>>& cs) {
    candidate_sets_.swap(cs);
    for (uint32_t i = 0; i < candidate_sets_.size(); ++i) {
      if (i == root_query_vertex_) continue;
      DLOG(INFO) << "query vertex " << i << ": " << candidate_sets_[i].size() << " candidates";
      ((TraverseOperator*)target_vertex_to_ops_[i])->setCandidateSets(&candidate_sets_[i]);
    }
  }

  inline const std::vector<VertexID>& getCandidateSet(QueryVertexID id) const {
    DCHECK_LT(id, candidate_sets_.size());
    return candidate_sets_[id];
  }

  inline OperatorTree cloneOperators() const { return operators_.clone(); }

  inline QueryVertexID getRootQueryVertexID() const { return root_query_vertex_; }

  inline bool isInCover(QueryVertexID id) const { return cover_table_[id] == 1; }

 private:
  inline void printOperatorChain(const Operator* root, const std::string& indent = "") const {
    LOG(INFO) << indent << root->toString();
    if (root->getNext() != nullptr) {
      printOperatorChain(root->getNext(), indent + "| ");
    }
  }

  TraverseOperator* newExpandEdgeOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                          const std::vector<int>& cover_table);
  TraverseOperator* newExpandKeyKeyVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex);
  TraverseOperator* newExpandSetVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex);
  TraverseOperator* newExpandSetToKeyVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex);
  TraverseOperator* newExpandIntoOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex);
  Operator* newOutputOperator();
};

}  // namespace circinus

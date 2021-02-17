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

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/operators.h"

namespace circinus {

class ExecutionPlan {
 private:
  const QueryGraph* query_graph_;
  std::vector<int> cover_table_;
  std::vector<QueryVertexID> matching_order_;

  std::vector<std::vector<VertexID>> candidate_sets_;

 public:
  ~ExecutionPlan() {
    for (auto op : operators_) {
      delete op;
    }
    operators_.clear();
  }

  void populatePhysicalPlan(const QueryGraph* g, const std::vector<QueryVertexID>& matching_order,
                            const std::vector<int>& cover_table);

  inline void setCandidateSets(std::vector<std::vector<VertexID>>& cs) {
    candidate_sets_.swap(cs);
    for (uint32_t i = 1; i < candidate_sets_.size(); ++i) {
      ((TraverseOperator*)target_vertex_to_ops_[i])->setCandidateSets(&candidate_sets_[i]);
    }
  }

 private:
  ExpandEdgeOperator* newExpandEdgeOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex);
  ExpandVertexOperator* newExpandVertexOperator(std::vector<QueryVertexID>& parents, QueryVertexID target_vertex);

  std::vector<Operator*> operators_;
  std::unordered_map<QueryVertexID, Operator*> target_vertex_to_ops_;
};

}  // namespace circinus

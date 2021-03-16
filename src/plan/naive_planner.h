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

#include <queue>
#include <tuple>
#include <vector>

#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/query_graph.h"
#include "plan/execution_plan.h"

namespace circinus {

class NaivePlanner {
 private:
  const QueryGraph* query_graph_;
  const std::vector<double>* candidate_cardinality_;

  WeightedBnB vertex_cover_solver_;
  ExecutionPlan plan_;

  std::vector<QueryVertexID> matching_order_;

 public:
  NaivePlanner(QueryGraph* query_graph, std::vector<double>* candidate_cardinality)
      : query_graph_(query_graph),
        candidate_cardinality_(candidate_cardinality),
        vertex_cover_solver_(query_graph, *candidate_cardinality) {}

  bool hasValidCandidate();

  const auto& getMatchingOrder() const { return matching_order_; }

  ExecutionPlan* generatePlan(const std::vector<QueryVertexID>& use_order = {});

  uint32_t analyzeDynamicCoreCoverEager(const std::vector<QueryVertexID>& use_order = {});
  std::tuple<uint32_t, uint32_t, uint32_t> analyzeDynamicCoreCoverMWVC();

 private:
  std::vector<QueryVertexID> generateMatchingOrder(const QueryGraph* g, const std::vector<int>& core_table,
                                                   QueryVertexID start_vertex);

  QueryVertexID selectStartingVertex(const std::vector<QueryVertexID>& cover);
};

}  // namespace circinus

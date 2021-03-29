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
#include <utility>
#include <vector>

#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/query_graph.h"
#include "plan/execution_plan.h"
#include "utils/profiler.h"

namespace circinus {

struct CoverNode {
  std::vector<QueryVertexID> cover_;
  uint64_t cover_bits_;
  uint32_t cardinality_;
  uint32_t parents_all_keys_cardinality_;
  std::vector<uint32_t> parents_;
};

class NaivePlanner {
 private:
  const QueryGraph* query_graph_;
  const std::vector<double>* candidate_cardinality_;
  WeightedBnB vertex_cover_solver_;
  ExecutionPlan plan_;

  std::vector<std::vector<CoverNode>> covers_;
  std::vector<std::vector<QueryVertexID>> to_intersect_vertices_;
  std::vector<QueryVertexID> matching_order_;
  std::vector<double> level_cost_;

 public:
  NaivePlanner(QueryGraph* query_graph, std::vector<double>* candidate_cardinality)
      : query_graph_(query_graph),
        candidate_cardinality_(candidate_cardinality),
        vertex_cover_solver_(query_graph, *candidate_cardinality) {}

  bool hasValidCandidate();

  const auto& getMatchingOrder() const { return matching_order_; }
  const auto& getCovers() const { return covers_; }
  const auto& getToIntersectVertices() const { return to_intersect_vertices_; }

  ExecutionPlan* generatePlan(const std::vector<QueryVertexID>& use_order = {}, Profiler* profiler = nullptr);
  ExecutionPlan* generateSamplePlan(const std::vector<QueryVertexID>& use_order = {}, Profiler* profiler = nullptr);
  ExecutionPlan* generatePlanWithEagerDynamicCover(const std::vector<QueryVertexID>& use_order = {},
                                               Profiler* profiler = nullptr);

ExecutionPlan* generatePlanWithSampleExecution (const std::vector<std::vector<double>>& cardinality, const std::vector<double>& level_cost, Profiler* profiler);
  void generateCoverNode (const std::vector<QueryVertexID>& use_order);

  std::pair<uint32_t, uint32_t> analyzeDynamicCoreCoverEager(const std::vector<QueryVertexID>& use_order = {});
  std::tuple<uint32_t, uint32_t, uint32_t> analyzeDynamicCoreCoverMWVC(
      const std::vector<QueryVertexID>& use_order = {});

 private:
  std::vector<QueryVertexID> generateMatchingOrder(const QueryGraph* g, const std::vector<int>& core_table,
                                                   QueryVertexID start_vertex);

  QueryVertexID selectStartingVertex(const std::vector<QueryVertexID>& cover);

  std::pair<uint32_t, uint32_t> analyzeDynamicCoreCoverEagerInner(const std::vector<int>& query_graph_cover);
  unordered_map<QueryVertexID, uint32_t> getDynamicCoreCoverEager(const std::vector<int>& query_graph_cover);
};

}  // namespace circinus

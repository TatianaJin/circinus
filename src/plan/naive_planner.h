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

#include <cmath>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/query_graph.h"
#include "plan/execution_plan.h"

namespace circinus {

struct CoverNode {
  std::vector<QueryVertexID> cover;
  uint64_t cover_bits;  // assume query size smaller than 64
  std::vector<uint32_t> parents;
};

class NaivePlanner {
 private:
  const QueryGraph* query_graph_;
  bool use_two_hop_traversal_;
  const std::vector<double>* candidate_cardinality_;
  ExecutionPlan plan_;

  std::vector<std::vector<CoverNode>> covers_;
  std::vector<QueryVertexID> matching_order_;
  std::vector<double> level_cost_;

  std::vector<double> logCardinality() {
    std::vector<double> log_cardinality((*candidate_cardinality_).size(), 0);
    for (uint32_t i = 0; i < (*candidate_cardinality_).size(); ++i) {
      log_cardinality[i] = log2((*candidate_cardinality_)[i]);
    }
    return log_cardinality;
  }

 public:
  NaivePlanner(QueryGraph* query_graph, std::vector<double>* candidate_cardinality)
      : query_graph_(query_graph), candidate_cardinality_(candidate_cardinality), plan_(GraphType::Normal) {}

  NaivePlanner(QueryGraph* query_graph, bool use_two_hop_traversal, std::vector<double>* candidate_cardinality,
               GraphType type)
      : query_graph_(query_graph),
        use_two_hop_traversal_(use_two_hop_traversal),
        candidate_cardinality_(candidate_cardinality),
        plan_(type) {}

  bool hasValidCandidate();

  void setCandidateCardinality(const std::vector<double>* candidates) { candidate_cardinality_ = candidates; }
  const auto& getMatchingOrder() const { return matching_order_; }
  const auto& getCovers() const { return covers_; }

  ExecutionPlan* generatePlan(const std::vector<QueryVertexID>& use_order = {});
  ExecutionPlan* generatePlanWithEagerDynamicCover(const std::vector<QueryVertexID>& use_order = {});
  ExecutionPlan* generatePlanWithoutCompression(const std::vector<QueryVertexID>& use_order = {});
  ExecutionPlan* generatePlanWithDynamicCover(const GraphBase* data_graph,
                                              const std::vector<CandidateSetView>* candidate_views);

  ExecutionPlan* generatePlanWithSampleExecution(const std::vector<std::vector<double>>& cardinality,
                                                 const std::vector<double>& level_cost);
  void generateCoverNode(const std::vector<std::vector<double>>& cardinality);
  void generateOrder(const std::vector<QueryVertexID>& use_order);

  std::pair<uint32_t, uint32_t> analyzeDynamicCoreCoverEager(const std::vector<QueryVertexID>& use_order = {});
  std::tuple<uint32_t, uint32_t, uint32_t> analyzeDynamicCoreCoverMWVC(
      const std::vector<QueryVertexID>& use_order = {});

 private:
  double estimateCardinality(const GraphBase* data_graph, const std::vector<CandidateSetView>* candidate_views,
                             uint64_t cover_bits, uint32_t level);

  void getCoverCC(QueryVertexID qid, QueryVertexID cc_id, const uint64_t cover_bits, std::vector<QueryVertexID>& cc);
  unordered_map<VertexID, double> dfsComputeCost(QueryVertexID qid, uint64_t cover_bits, std::vector<bool>& visited,
                                                 const GraphBase* data_graph,
                                                 const std::vector<CandidateSetView>* candidate_views,
                                                 const std::vector<QueryVertexID>& cc, bool with_traversal = false);
  std::vector<QueryVertexID> generateMatchingOrder(const QueryGraph* g, const std::vector<int>& core_table,
                                                   QueryVertexID start_vertex);

  QueryVertexID selectStartingVertex(const std::vector<QueryVertexID>& cover);

  std::pair<uint32_t, uint32_t> analyzeDynamicCoreCoverEagerInner(const std::vector<int>& query_graph_cover);
  unordered_map<QueryVertexID, uint32_t> getDynamicCoreCoverEager(const std::vector<int>& query_graph_cover);

  std::vector<int> getCoverByOrder() const;
};

}  // namespace circinus

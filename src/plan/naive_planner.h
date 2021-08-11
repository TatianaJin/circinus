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
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/query_graph.h"
#include "plan/execution_plan.h"

namespace circinus {

/**
 * The following functions should be overrode to implement a different compression strategy:
 *  - For pruning cover search space
 *    - generateAnchorCovers: compute the anchor covers from which sequences of covers are generated to populate the
 *                            search space
 *  - For estimating the cardinaltiy
 *    - getMinimalConectedSubgraphWithAllKeys: compute the connected subgraph on which cardinality is estimated
 *    - dfsComputeCost: compute the cost using a dfs tree on the connected subgraph
 *    - estimateCardinality: estimate the cardinality of cover embeddings of a subquery
 *  - For estimating the computation costs
 *    - estimateExpandCost: estimate the computation cost based on the cardinality estimation
 * The following functions should be overrode to implement a different parallelization strategy:
 *  - selectParallelizingQueryVertex
 *  - getParallelizingQueryVertexWeights
 */
class NaivePlanner {
 private:
  struct CoverNode {
    std::vector<QueryVertexID> cover;
    uint64_t cover_bits = 0;  // assume query size smaller than 64
    std::vector<uint32_t> parents;

    inline std::vector<int> getCoverTable(QueryVertexID size) const {
      std::vector<int> select_cover(size, 0);
      for (uint32_t j = 0; j < size; ++j) {
        if (cover_bits >> j & 1) {
          select_cover[j] = 1;
        }
      }
      return select_cover;
    }

    inline std::string getCoverTableString(QueryVertexID size) const {
      std::stringstream ss;
      ss << (cover_bits & 1);
      for (uint32_t j = 1; j < size; ++j) {
        ss << ',' << (cover_bits >> j & 1);
      }
      return ss.str();
    }
  };

  const QueryGraph* query_graph_;
  bool use_two_hop_traversal_ = false;
  const std::vector<double> candidate_cardinality_;
  ExecutionPlan plan_;

  std::vector<std::vector<CoverNode>> covers_;  // covers for each level, the search space for compression
  std::vector<QueryVertexID> matching_order_;
  std::vector<double> level_cost_;

 public:
  NaivePlanner(QueryGraph* query_graph, std::vector<double>&& candidate_cardinality)
      : query_graph_(query_graph), candidate_cardinality_(std::move(candidate_cardinality)), plan_(GraphType::Normal) {}

  NaivePlanner(QueryGraph* query_graph, bool use_two_hop_traversal, std::vector<double>&& candidate_cardinality,
               GraphType type)
      : query_graph_(query_graph),
        use_two_hop_traversal_(use_two_hop_traversal),
        candidate_cardinality_(std::move(candidate_cardinality)),
        plan_(type) {}

  const auto& getMatchingOrder() const { return matching_order_; }
  const auto& getCovers() const { return covers_; }

  // FIXME(tatiana): tidy-up the interface and implementations
  ExecutionPlan* generatePlan(const std::vector<QueryVertexID>& use_order = {});
  ExecutionPlan* generatePlanWithEagerDynamicCover(const std::vector<QueryVertexID>& use_order = {});
  ExecutionPlan* generatePlanWithoutCompression(const std::vector<QueryVertexID>& use_order = {});
  /**
   * If candidate_views is nullptr, estimate the cardinality of compressed groups by product of candidate cardinality;
   * otherwise, do a dfs traversal based on the candidate view to compute a tighter estimation.
   */
  ExecutionPlan* generatePlanWithDynamicCover(const GraphBase* data_graph,
                                              const std::vector<CandidateSetView>* candidate_views);

  ExecutionPlan* generatePlanWithSampleExecution(const std::vector<std::vector<double>>& cardinality,
                                                 const std::vector<double>& level_cost);
  /** Generates the pruned search space of dynamic cover sequences */
  void generateCoverNode(const std::vector<std::vector<double>>& cardinality);
  const std::vector<QueryVertexID>& generateOrder(const GraphBase* data_graph, const GraphMetadata& metadata,
                                                  const std::vector<CandidateSetView>* candidate_views,
                                                  const std::vector<VertexID>& candidate_cardinality,
                                                  OrderStrategy order_strategy,
                                                  const std::vector<QueryVertexID>& use_order = {});

  inline std::pair<QueryVertexID, std::vector<double>> computeParallelVertexWeights(
      const GraphBase* data_graph, const std::vector<CandidateSetView>* candidate_views,
      const std::vector<QueryVertexID>& parallelizing_qv_candidates) {
    auto cover_bits = plan_.getQueryCoverBits();
    DCHECK_NE(cover_bits, 0);
    auto parallelizing_qv = selectParallelizingQueryVertex(cover_bits, parallelizing_qv_candidates);
    CHECK_NE(parallelizing_qv, DUMMY_QUERY_VERTEX);
    return std::make_pair(parallelizing_qv, getParallelizingQueryVertexWeights(parallelizing_qv, data_graph,
                                                                               candidate_views, cover_bits));
  }

  std::pair<uint32_t, uint32_t> analyzeDynamicCoreCoverEager(const std::vector<QueryVertexID>& use_order = {});
  std::tuple<uint32_t, uint32_t, uint32_t> analyzeDynamicCoreCoverMWVC(
      const std::vector<QueryVertexID>& use_order = {});

 private:
  /* Start of implementations of compression strategy */

  std::vector<std::vector<int>> generateAnchorCovers(const std::vector<QueryVertexID>& subquery_vertices,
                                                     const std::vector<double>& cardinality);
  /** Alternative for generateAnchorCovers */
  std::vector<std::vector<int>> generateAnchorCoversSettingParentAsKey(
      const std::vector<QueryVertexID>& subquery_vertices, const std::vector<double>& cardinality,
      const std::vector<uint32_t>& order_index,
      const std::vector<std::vector<QueryVertexID>>& parent_vertices_per_level,
      const std::vector<unordered_set<QueryVertexID>>& existing_vertices);

  /** Computes a connected subgraph of cover vertices given in the cover bits for the subquery at the given level.
   *
   * First find connected components (CCs) in the cover of the subquery at the given level, and then greedily pick
   * non-cover vertices with smallest cardinality to connect CCs into one connected subgraph.
   */
  void getMinimalConectedSubgraphWithAllKeys(std::vector<QueryVertexID>& cc, uint32_t level, uint64_t cover_bits);

  /** Estimates the cardinality of compressed groups for the subquery given the cover.
   *
   * First computes a minimal connected subgraph containing all cover vertices, and estimate the cardinality by path
   * weights based on a dfs tree.
   */
  double estimateCardinality(const GraphBase* data_graph, const std::vector<CandidateSetView>* candidate_views,
                             uint64_t cover_bits, uint32_t level);

  inline double estimateCardinality(const GraphBase* data_graph, const std::vector<CandidateSetView>* candidate_views,
                                    const CoverNode& cover_node, uint32_t level) {
    if (candidate_views == nullptr) {
      double car = 1;
      for (QueryVertexID qid : cover_node.cover) {
        car *= candidate_cardinality_[qid];
      }
      return car;
    }
    return estimateCardinality(data_graph, candidate_views, cover_node.cover_bits, level);
  }

  double estimateExpandCost(const GraphBase* data_graph, const std::vector<CandidateSetView>* candidate_views,
                            const std::vector<std::vector<double>>& car,
                            const unordered_set<QueryVertexID>& existing_vertices, uint32_t level, uint32_t idx,
                            uint32_t parent);

  /**
   * We estimate the cardinality of cover embeddings by first-order traversal (virtually) on the bipartite graphs of
   * candidate sets. The idea is similar to that of computing path weights in CFL.
   */
  template <typename Hashmap = unordered_map<VertexID, double>>
  unordered_map<VertexID, double> dfsComputeCost(QueryVertexID qid, uint64_t cover_bits, std::vector<bool>& visited,
                                                 std::vector<bool>& visited_cc, const GraphBase* data_graph,
                                                 const std::vector<CandidateSetView>* candidate_views,
                                                 const std::vector<QueryVertexID>& cc, bool with_traversal = false);

  /* End of implementations of compression strategy */

  /* start of parallelization strategy */

  QueryVertexID selectParallelizingQueryVertex(uint64_t cover_bits,
                                               const std::vector<QueryVertexID>& parallelizing_qv_candidates);

  std::vector<double> getParallelizingQueryVertexWeights(QueryVertexID partition_qv, const GraphBase* data_graph,
                                                         const std::vector<CandidateSetView>* candidate_views,
                                                         uint64_t cover_bits);
  /* end of parallelization strategy */

  bool hasValidCandidate();

  inline std::vector<double> logCardinality() {
    std::vector<double> log_cardinality(candidate_cardinality_.size(), 0);
    for (uint32_t i = 0; i < candidate_cardinality_.size(); ++i) {
      log_cardinality[i] = log2(candidate_cardinality_[i]);
    }
    return log_cardinality;
  }

  inline bool addCover(CoverNode& new_cover_node, uint32_t level) {
    {  // debug log
      std::stringstream ss;
      for (auto x : new_cover_node.cover) {
        ss << x << ", ";
      }
      DLOG(INFO) << level << " " << ss.str();
    }
    for (const auto& cover_node : covers_[level]) {
      if (cover_node.cover_bits == new_cover_node.cover_bits) {
        return true;
      }
    }
    covers_[level].emplace_back(new_cover_node);
    return false;
  }

  void getCoverCC(QueryVertexID qid, QueryVertexID cc_id, const uint64_t cover_bits, std::vector<QueryVertexID>& cc);

  std::vector<QueryVertexID> generateMatchingOrder(const QueryGraph* g, const std::vector<int>& core_table,
                                                   QueryVertexID start_vertex);

  // FIXME(tatiana): reuse the implementation of ordering strategy
  QueryVertexID selectStartingVertex(const std::vector<QueryVertexID>& cover);

  std::pair<uint32_t, uint32_t> analyzeDynamicCoreCoverEagerInner(const std::vector<int>& query_graph_cover);
  unordered_map<QueryVertexID, uint32_t> getDynamicCoreCoverEager(const std::vector<int>& query_graph_cover);

  std::vector<int> getCoverByOrder() const;
};

}  // namespace circinus

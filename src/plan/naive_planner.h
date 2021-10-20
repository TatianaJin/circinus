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

#include <algorithm>
#include <cmath>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/partial_order.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "plan/execution_plan.h"
#include "utils/hashmap.h"

namespace circinus {

/**
 * The following functions should be overrode to implement a different compression strategy:
 *  - For pruning cover search space
 *    - generateAnchorCovers: compute the anchor covers from which sequences of covers
 *                            are generated to populate the search space
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
  /** The search space of the dynamic cover sequences is a set of DAGs, in which CoverNode is a node.
   *
   * A parent of this is a CoverNode that corresponds to the next-level subquery, and whose cover has at most one more
   * vertex and is compatible with the this cover.
   */
  struct CoverNode {
    uint64_t cover_bits = 0;  // assume query size smaller than 64
    std::vector<uint32_t> parents;

    inline std::vector<int> getCoverTable(QueryVertexID size) const {
      std::vector<int> select_cover(size, 0);
      for (uint32_t j = 0; j < size; ++j) {
        select_cover[j] = (cover_bits >> j & 1);
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

  const QueryGraph* const query_graph_;
  QueryGraph ordered_query_graph_;
  bool use_two_hop_traversal_ = false;
  std::vector<double> candidate_cardinality_;
  ExecutionPlan plan_;

  std::vector<std::vector<CoverNode>> covers_;  // covers for each level, the search space for compression
  std::vector<QueryVertexID> matching_order_;

  std::unique_ptr<VertexRelationship> vertex_relationship_ = nullptr;

 public:
  NaivePlanner(QueryGraph* query_graph, bool use_two_hop_traversal, std::vector<double>&& candidate_cardinality,
               GraphType type)
      : query_graph_(query_graph),
        use_two_hop_traversal_(use_two_hop_traversal),
        candidate_cardinality_(std::move(candidate_cardinality)),
        plan_(type) {
    CHECK_LE(query_graph->getNumVertices(), 64);
  }

  NaivePlanner(QueryGraph* query_graph, GraphType type) : query_graph_(query_graph), plan_(type) {}

  inline void setVertexEquivalence(const VertexEquivalence& ve) {
    vertex_relationship_ = std::make_unique<VertexRelationship>(ve);
    plan_.setVertexRelationship(vertex_relationship_.get());
  }

  /** Generates a plan with a static cover for partial match compression. */
  ExecutionPlan* generatePlan(const PartialOrder* = nullptr);

  /** Generates a plan with no compression. */
  ExecutionPlan* generatePlanWithoutCompression();

  /**
   * If candidate_views is nullptr, estimate the cardinality of compressed groups by product of candidate cardinality;
   * otherwise, do a dfs traversal based on the candidate view to compute a tighter estimation.
   */
  ExecutionPlan* generatePlanWithDynamicCover(const GraphBase* data_graph,
                                              const std::vector<CandidateSetView>* candidate_views,
                                              const PartialOrder* = nullptr);

  ExecutionPlan* generatePlanWithSampleExecution(const std::vector<std::vector<double>>& cardinality,
                                                 const std::vector<double>& level_cost);

  /** Generates query vertex matching order.
   *
   * Candidate views are required for CFL, DAF, TSO and GQL.
   * TODO(tatiana): handle case when candidate views are not feasible to compute.
   */
  const std::vector<QueryVertexID>& generateOrder(const GraphBase* data_graph, const GraphMetadata& metadata,
                                                  const std::vector<CandidateSetView>* candidate_views,
                                                  const std::vector<VertexID>& candidate_cardinality,
                                                  OrderStrategy order_strategy, QueryVertexID seed_qv,
                                                  const std::vector<QueryVertexID>* use_order = nullptr);

  /** Generates query vertex matching order for online query.
   *
   * If seed_qv is a query vertex, generate an order starting from the seed vertex.
   */
  const std::vector<QueryVertexID>& generateOrder(QueryVertexID seed_qv, OrderStrategy os, const PartialOrder* po);

  /** Select a query vertex to parallelize the plan and compute the weights (estimation of the associated backtracking
   * search space) for all its candidates. */
  inline std::pair<QueryVertexID, std::vector<double>> computeParallelVertexWeights(
      const GraphBase* data_graph, const std::vector<CandidateSetView>& candidate_views,
      const std::vector<QueryVertexID>& parallelizing_qv_candidates) {
    auto cover_bits = plan_.getQueryCoverBits();
    CHECK_NE(cover_bits, 0);
    auto parallelizing_qv = selectParallelizingQueryVertex(cover_bits, parallelizing_qv_candidates);
    CHECK_NE(parallelizing_qv, DUMMY_QUERY_VERTEX);
    return std::make_pair(parallelizing_qv, getParallelizingQueryVertexWeights(parallelizing_qv, data_graph,
                                                                               candidate_views, cover_bits));
  }

  /* start of parallelization strategy */

  std::vector<double> getParallelizingQueryVertexWeights(QueryVertexID partition_qv, const GraphBase* data_graph,
                                                         const std::vector<CandidateSetView>& candidate_views,
                                                         uint64_t cover_bits);

  QueryVertexID selectParallelizingQueryVertex(uint64_t cover_bits,
                                               const std::vector<QueryVertexID>& parallelizing_qv_candidates,
                                               double parallelization_threshold = 5e5);

  /* end of parallelization strategy */

  inline const std::vector<QueryVertexID>& getMatchingOrder() const { return matching_order_; }

 private:
  /* Start of implementations of compression strategy */

  /** Generates the pruned search space of dynamic cover sequences */
  void generateCoverNode(const std::vector<std::vector<double>>& cardinality);

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
  double estimateCardinality(const GraphBase* data_graph, const std::vector<CandidateSetView>& candidate_views,
                             uint64_t cover_bits, uint32_t level);

  inline double estimateCardinality(const std::vector<double>& cardinality, const GraphBase* data_graph,
                                    const std::vector<CandidateSetView>* candidate_views, const CoverNode& cover_node,
                                    uint32_t level) {
    if (candidate_views == nullptr) {
      return estimateCardinalityByProduct(cardinality, cover_node.cover_bits, level);
    }
    return estimateCardinality(data_graph, *candidate_views, cover_node.cover_bits, level);
  }

  inline double estimateCardinalityByProduct(const std::vector<double>& cardinality, uint64_t cover_bits,
                                             uint32_t level) const {
    double car = 1;
    for (uint32_t i = 0; i <= level; ++i) {
      auto qid = matching_order_[i];
      if (cover_bits >> qid & 1) {
        car *= cardinality[qid];
      }
    }
    return car;
  }

  double estimateExpandCost(const GraphBase* data_graph, const std::vector<CandidateSetView>* candidate_views,
                            std::vector<unordered_map<uint64_t, double>>& car,
                            const unordered_set<QueryVertexID>& existing_vertices, uint32_t level, uint32_t idx,
                            uint32_t parent);

  /**
   * We estimate the cardinality of cover embeddings by first-order traversal (virtually) on the bipartite graphs of
   * candidate sets. The idea is similar to that of computing path weights in CFL.
   */
  unordered_map<VertexID, double> dfsComputeCost(QueryVertexID qid, uint64_t cover_bits, std::vector<bool>& visited,
                                                 std::vector<bool>& visited_cc, const GraphBase* data_graph,
                                                 const std::vector<CandidateSetView>& candidate_views,
                                                 const std::vector<QueryVertexID>& cc, bool with_traversal = false);

  double intersectionCountCoeff(const std::vector<double>& cardinality,
                                const unordered_set<QueryVertexID>& existing_vertices, uint32_t level, uint32_t idx,
                                uint32_t parent);

  /* End of implementations of compression strategy */

  bool hasValidCandidate() const;

  inline std::vector<double> logCardinality() {
    std::vector<double> log_cardinality(candidate_cardinality_.size(), 0);
    for (uint32_t i = 0; i < candidate_cardinality_.size(); ++i) {
      log_cardinality[i] = log2(candidate_cardinality_[i]);
    }
    return log_cardinality;
  }

  void setVertexWeightByDegreeConstraints(std::vector<double>& weights, const PartialOrder* po) const;

  inline bool addCover(CoverNode& new_cover_node, uint32_t level) {
    DLOG(INFO) << level << " " << new_cover_node.getCoverTableString(matching_order_.size()) << " target "
               << matching_order_[level];
    for (const auto& cover_node : covers_[level]) {
      if (cover_node.cover_bits == new_cover_node.cover_bits) {
        return true;
      }
    }
    covers_[level].emplace_back(new_cover_node);
    return false;
  }

  void getCoverCC(QueryVertexID qid, QueryVertexID cc_id, const uint64_t cover_bits,
                  std::vector<QueryVertexID>& cc) const;

  std::pair<unordered_map<QueryVertexID, uint32_t>, std::vector<double>> traceBackCoverPath(
      const std::vector<std::vector<double>>& costs_car, const std::vector<std::vector<uint32_t>>& pre, int best_idx);

  void logCoverSpace();

  inline double getCardinality(std::vector<unordered_map<uint64_t, double>>& car, uint32_t level, uint64_t cover_bits,
                               const GraphBase* data_graph, const std::vector<CandidateSetView>* candidate_views) {
    if
      constexpr(false) {
        // share costs among identical covers across levels
        auto pos = car.front().find(cover_bits);
        if (pos == car.front().end()) {
          auto cardinality = estimateCardinality(data_graph, *candidate_views, cover_bits, level);
          pos = car.front().insert({cover_bits, cardinality}).first;
        }
        return pos->second;
      }
    // seperate cover cost for each level
    auto pos = car[level].find(cover_bits);
    if (pos == car[level].end()) {
      auto cardinality = estimateCardinality(data_graph, *candidate_views, cover_bits, level);
      if
        constexpr(false && level > 1) {
          auto last_pos = car[level - 1].find(cover_bits);
          if (last_pos != car[level - 1].end()) {
            cardinality = std::min(cardinality, last_pos->second);
          }
        }
      pos = car[level].insert({cover_bits, cardinality}).first;
    }
    return pos->second;
  }
};

}  // namespace circinus

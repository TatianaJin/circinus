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

#include "plan/planner.h"

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "algorithms/centrality.h"
#include "ops/logical/compressed_input.h"
#include "plan/backtracking_plan.h"
#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/query_utils.h"

namespace circinus {

/**
 * For now, in the case of partitioned graph, all graph partitions share the same candidate pruning plan.
 */
CandidatePruningPlan* Planner::generateCandidatePruningPlan() {
  auto& q = query_context_->query_graph;
  auto& strategy = query_context_->query_config.candidate_pruning_strategy;
  switch (strategy) {
  case CandidatePruningStrategy::None: {
    candidate_pruning_plan_.setFinished();
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::Adaptive: {
    LOG(FATAL) << "Not implemented yet";
    // TODO(tatiana): adapt to data graph distribution
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::LDF: {
    LOG(INFO) << "Prune candidate by LDF";
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.setPartitionResult(toPartitionCandidates());
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::NLF: {
    LOG(INFO) << "Prune candidate by NLF";
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
    candidate_pruning_plan_.setPartitionResult(toPartitionCandidates());
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::DAF: {
    candidate_pruning_plan_.newLDFScan(q);
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::CFL: {
    // auto& metadata = *query_context_->graph_metadata;
    // auto order = std::make_unique<CFLOrder>();
    // candidate_pruning_plan_.newLDFScan(q, order->getTopThree(metadata, &q));
    // order_ = std::move(order);
    candidate_pruning_plan_.newLDFScan(q);  // now we generate the candidates for all query vertices at once
    candidate_pruning_plan_.newNLFFilter(q, true);
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::TSO: {
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::GQL: {
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
    return &candidate_pruning_plan_;
  }
  }
  return &candidate_pruning_plan_;
}

// TODO(engineering): classes in plan should not know classes in exec
CandidatePruningPlan* Planner::updateCandidatePruningPlan(const CandidateResult* result) {
  CHECK(result != nullptr);
  auto part_cardinality = result->getCandidateCardinality();
  auto& q = query_context_->query_graph;
  auto& strategy = query_context_->query_config.candidate_pruning_strategy;
  auto n_qvs = part_cardinality[0].size();
  std::vector<VertexID> cardinality(n_qvs, 0);
  for (uint32_t i = 0; i < part_cardinality.size(); ++i) {
    for (uint32_t j = 0; j < n_qvs; ++j) {
      cardinality[j] += part_cardinality[i][j];
    }
  }
  if (strategy == CandidatePruningStrategy::LDF || strategy == CandidatePruningStrategy::NLF) {
    for (uint32_t i = 0; i < cardinality.size(); ++i) {
      DLOG(INFO) << " |C(v" << i << ")|: " << cardinality[i];
    }
    candidate_pruning_plan_.setFinished();
    return &candidate_pruning_plan_;
  }
  auto phase = candidate_pruning_plan_.completePhase();
  if (phase == 2) {
    // TODO(tatiana): now we skip phase 2 as it is not easy to parallelize forward construction due to set union
    phase = candidate_pruning_plan_.completePhase();
  }
  if (phase == 3) {
    auto& metadata = *query_context_->graph_metadata;
    switch (strategy) {
    case CandidatePruningStrategy::DAF: {
      candidate_pruning_plan_.newDAFFilter(&q, metadata, cardinality);
      break;
    }
    case CandidatePruningStrategy::CFL: {
      candidate_pruning_plan_.newCFLFilter(&q, metadata, cardinality);
      break;
    }
    case CandidatePruningStrategy::GQL: {
      candidate_pruning_plan_.newGQLFilter(&q);
      break;
    }
    case CandidatePruningStrategy::TSO: {
      candidate_pruning_plan_.newTSOFilter(&q, metadata, cardinality);
      break;
    }
    default:
      candidate_pruning_plan_.setFinished();
    }
    return &candidate_pruning_plan_;
  }
  candidate_pruning_plan_.setPartitionResult(toPartitionCandidates());
  for (uint32_t i = 0; i < cardinality.size(); ++i) {
    DLOG(INFO) << " |C(v" << i << ")|: " << cardinality[i];
  }
  candidate_pruning_plan_.setFinished();
  return &candidate_pruning_plan_;
}

std::vector<std::vector<VertexID>> Planner::estimateCardinality() const {
  auto metadata = query_context_->graph_metadata;
  std::vector<std::vector<VertexID>> ret;
  ret.reserve(metadata->numPartitions());
  if (toPartitionCandidates()) {
    for (uint32_t i = 0; i < metadata->numPartitions(); ++i) {
      ret.emplace_back(estimateCardinalityInner(&metadata->getPartition(i)));
    }
  } else {
    ret.emplace_back(estimateCardinalityInner(metadata));
  }
  return ret;
}

std::vector<QueryVertexID> Planner::getPartitioningQueryVertices() {
  switch (query_context_->query_config.pqv_strategy) {
  case PQVStrategy::None: {
    return {};
  }
  case PQVStrategy::ClosenessCentrality: {
    return getPartitioningQueryVerticesByClosenessCentrality();
  }
  }
  return {};
}

std::vector<QueryVertexID> Planner::getPartitioningQueryVerticesByCover() {
  auto& q = query_context_->query_graph;

  // we favor the query vertices that are used as the compression key for task data partitioning
  std::vector<uint32_t> qv_frequency_in_cover(q.getNumVertices(), 0);
  for (auto plan : backtracking_plan_->getPlans()) {
    for (QueryVertexID i = 0; i < q.getNumVertices(); ++i) {
      qv_frequency_in_cover[i] += plan->isInCover(i);
    }
  }
  std::vector<QueryVertexID> order_by_cover_occurrence(q.getNumVertices());
  std::iota(order_by_cover_occurrence.begin(), order_by_cover_occurrence.end(), 0);
  std::sort(order_by_cover_occurrence.begin(), order_by_cover_occurrence.end(),
            [&qv_frequency_in_cover](QueryVertexID a, QueryVertexID b) {
              return qv_frequency_in_cover[a] > qv_frequency_in_cover[b];
            });

  uint32_t common_cover_vertex_count = 0;
  std::stringstream ss;
  for (auto v : order_by_cover_occurrence) {
    if (qv_frequency_in_cover[v] == backtracking_plan_->getPlans().size()) {
      ++common_cover_vertex_count;
      ss << ' ' << v;
    } else {
      break;
    }
  }
  LOG(INFO) << "Common cover vertex: [" << ss.str() << ']';
  if (common_cover_vertex_count != 0) {
    // TODO(tatiana): pick the indicator by considering the closeness centrality?
  }

  // TODO(tatiana): now we use only one partitioning query vertex to test the pipeline, later may consider
  // parallelism for deciding partitioning vertices
  order_by_cover_occurrence.resize(1);

  return order_by_cover_occurrence;
}

std::vector<QueryVertexID> Planner::getPartitioningQueryVerticesByClosenessCentrality() {
  auto& q = query_context_->query_graph;
  auto centrality_gen = ClosenessCentrality(&q);
  auto centrality = centrality_gen.getClosenessCentrality();

  std::vector<QueryVertexID> order_by_closeness_centrality(q.getNumVertices());
  std::iota(order_by_closeness_centrality.begin(), order_by_closeness_centrality.end(), 0);
  std::sort(order_by_closeness_centrality.begin(), order_by_closeness_centrality.end(),
            [&centrality](QueryVertexID a, QueryVertexID b) { return centrality[a] > centrality[b]; });

  if (verbosePlannerLog()) {
    std::stringstream ss;
    for (auto v : order_by_closeness_centrality) {
      ss << ' ' << v << ": " << centrality[v] << ",  ";
    }
    LOG(INFO) << "Query vertex centrality: [" << ss.str() << ']';
  }

  order_by_closeness_centrality.resize(1);

  return order_by_closeness_centrality;
}

void Planner::exhaustivePartitionPlan(std::vector<std::vector<CandidateScope>>& partitioned_plans, uint32_t level,
                                      uint32_t partition_num, const std::vector<QueryVertexID>& partitioning_qv,
                                      std::vector<CandidateScope> partitioned_plan) {
  if (level == partitioning_qv.size()) {
    if (verbosePlannerLog()) {  // debug log
      std::stringstream ss;
      for (auto x : partitioned_plan) {
        x.print(ss);
      }
      LOG(INFO) << ss.str();
    }
    partitioned_plans.emplace_back(std::move(partitioned_plan));
    return;
  }

  for (uint32_t i = 0; i < partition_num; ++i) {
    std::vector<CandidateScope> new_partitioned_plan(partitioned_plan.begin(), partitioned_plan.end());
    new_partitioned_plan[partitioning_qv[level]].usePartition(i);
    exhaustivePartitionPlan(partitioned_plans, level + 1, partition_num, partitioning_qv, new_partitioned_plan);
  }
}

std::vector<std::vector<CandidateScope>> Planner::generatePartitionedScopes(
    uint32_t partition_num, const std::vector<QueryVertexID>& partitioning_qvs,
    PartitionedCandidateResult* partitioned_candidates) {
  std::vector<std::vector<CandidateScope>> partitioned_plans;
  uint32_t n_qvs = query_context_->query_graph.getNumVertices();
  exhaustivePartitionPlan(partitioned_plans, 0, partition_num, partitioning_qvs, std::vector<CandidateScope>(n_qvs));
  return partitioned_plans;
}

// TODO(tatiana): remove obsolete functions
std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> Planner::generatePartitionedPlans(
    const std::vector<QueryVertexID>& partitioning_qv) {
  // CHECK_EQ(partitioning_qv.size(), 1);  // consider the locality indicator only, no. tasks equal to no. partitions
  DLOG(INFO) << "partitioning_qv = " << partitioning_qv.front();
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> ret;
  ret.reserve(backtracking_plan_->getPlans().size() * ((1 << partitioning_qv.size()) - 1));

  auto n_qvs = query_context_->query_graph.getNumVertices();

  for (uint32_t i = 0; i < backtracking_plan_->getPlans().size(); ++i) {
    ret.emplace_back(i, std::vector<CandidateScope>(n_qvs));
    ret.back().second[partitioning_qv.front()].usePartition(backtracking_plan_->getPlan(i)->getPartitionId());
  }
  LOG(INFO) << "Partition qeury vertex size: " << partitioning_qv.size() << ",  partition plan size:" << ret.size();
  return ret;
}

void Planner::newInputOperators(ExecutionPlan* logical_plan, const std::vector<QueryVertexID>* partitioning_qvs) {
  if (partitioning_qvs == nullptr || partitioning_qvs->empty()) {
    backtracking_plan_->addInputOperator(std::make_unique<LogicalCompressedInputOperator>(
        logical_plan->getMatchingOrder().front(), logical_plan->inputAreKeys()));
  } else {
    backtracking_plan_->addInputOperator(std::make_unique<PartitionedLogicalCompressedInputOperator>(
        &query_context_->query_graph, logical_plan->inputAreKeys(), logical_plan->getMatchingOrder(),
        *partitioning_qvs));
  }
}

std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> Planner::generateLogicalPlans(
    const std::vector<QueryVertexID>& partitioning_qvs, const std::vector<std::vector<CandidateScope>>& partitions,
    PartitionedCandidateResult* result) {
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> partition_plans;
  partition_plans.reserve(partitions.size());
  bool intra_plan_using_primary_pqv = query_context_->query_config.intra_partition_plan;
  if (partitioning_qvs.empty()) intra_plan_using_primary_pqv = false;
  if (intra_plan_using_primary_pqv) {  // generate a plan for each different scope on primary pqv
    auto primary_pqv = partitioning_qvs.front();
    uint32_t n_qvs = query_context_->query_graph.getNumVertices();
    auto candidate_cardinality = result->getCandidateCardinality();
    unordered_map<uint32_t, uint32_t> scope_partition_to_plan_idx;
    int global_plan_idx = -1;  // in case a global plan is needed
    for (auto& scopes : partitions) {
      // assume the candidates of the primary pqv are split for each graph partition
      CHECK(scopes[primary_pqv].getType() == CandidateScopeType::Partition) << (uint32_t)scopes[primary_pqv].getType();
      auto graph_partition = scopes[primary_pqv].getPartition();
      auto plan_pos = scope_partition_to_plan_idx.find(graph_partition);
      if (plan_pos == scope_partition_to_plan_idx.end()) {  // plan not generated yet, generate intra partition plan
        LOG(INFO) << "Generate Plan " << backtracking_plan_->getPlans().size() << " for intra-partition "
                  << graph_partition;
        std::vector<CandidateScope> plan_scopes(n_qvs, scopes[primary_pqv]);
        // std::vector<CandidateScope> plan_scopes(n_qvs);
        // plan_scopes[primary_pqv] = scopes[primary_pqv];
        auto& stats = candidate_cardinality[graph_partition];
        auto candidate_views = result->getCandidatesByScopes(plan_scopes);
        auto plan = generateLogicalExecutionPlan(stats, primary_pqv, &candidate_views, &partitioning_qvs);
        if (plan == nullptr) {
          if (planners_.size() == backtracking_plan_->getPlans().size() + 1) {  // need a global plan
            planners_.pop_back();
            if (global_plan_idx == -1) {
              auto candidate_views = result->getCandidates();
              plan = generateLogicalExecutionPlan(getCardinality(candidate_views), DUMMY_QUERY_VERTEX, &candidate_views,
                                                  &partitioning_qvs);
              CHECK(plan != nullptr);  // the entire problem does not have a match if nullptr
              global_plan_idx = backtracking_plan_->addPlan(plan);
            }
            plan_pos = scope_partition_to_plan_idx.insert({graph_partition, global_plan_idx}).first;
          } else {
            continue;
          }
        } else {
          plan_pos = scope_partition_to_plan_idx.insert({graph_partition, backtracking_plan_->addPlan(plan)}).first;
        }
      }
      // assign plan to scope
      partition_plans.emplace_back(plan_pos->second, scopes);
    }
    return partition_plans;
  }
  // generate a different plan for each backtracking scope
  for (auto& scopes : partitions) {
    LOG(INFO) << "Generate Plan " << backtracking_plan_->getPlans().size();
    auto candidate_views = result->getCandidatesByScopes(scopes);
    auto stats = getCardinality(candidate_views);
    auto plan = generateLogicalExecutionPlan(stats, DUMMY_QUERY_VERTEX, &candidate_views, &partitioning_qvs);
    if (plan == nullptr) {
      planners_.resize(backtracking_plan_->getPlans().size());
      continue;
    }
    partition_plans.emplace_back(backtracking_plan_->addPlan(plan), scopes);
  }
  return partition_plans;
}

void Planner::parallelizePartitionedPlans(
    const std::vector<QueryVertexID>& partitioning_qvs,
    std::vector<std::pair<uint32_t, std::vector<CandidateScope>>>* partitioned_plans,
    PartitionedCandidateResult* result) {
  // TODO(tatiana): handle when candidate views are not available
  /* get weights for each plan partition */
  double max_weight = 0;
  std::vector<std::pair<QueryVertexID, std::vector<double>>> plan_weights;
  plan_weights.reserve(partitioned_plans->size());
  for (auto& partition : *partitioned_plans) {
    auto candidate_views = result->getCandidatesByScopes(partition.second);
    auto[parallel_qv, weights] = planners_[partition.first]->computeParallelVertexWeights(
        query_context_->data_graph, candidate_views, partitioning_qvs);
    std::stringstream ss;
    bool parallel_qv_is_pqv = false;
    for (auto v : partitioning_qvs) {
      ss << ' ' << v;
      if (v == parallel_qv) parallel_qv_is_pqv = true;
    }
    if (parallel_qv_is_pqv || backtracking_plan_->getPlan(partition.first)->getMatchingOrder()[0] == parallel_qv) {
      if (verbosePlannerLog()) {
        LOG(INFO) << "parallel_qv = " << parallel_qv << " partitioning_qvs [" << ss.str() << " ]";
      }
    } else {  // replace existing logical input operator
      std::vector<QueryVertexID> pqvs(partitioning_qvs.size() + 1);
      std::copy(partitioning_qvs.begin(), partitioning_qvs.end(), pqvs.begin());
      pqvs.back() = parallel_qv;
      auto logical_plan = backtracking_plan_->getPlan(partition.first);
      backtracking_plan_->replaceInputOperator(
          partition.first,
          std::make_unique<PartitionedLogicalCompressedInputOperator>(
              &query_context_->query_graph, logical_plan->inputAreKeys(), logical_plan->getMatchingOrder(), pqvs));
      LOG(INFO) << "Regenerate input operator as parallel qv is not pqv. parallel_qv = " << parallel_qv
                << " partitioning_qvs [" << ss.str() << " ]";
    }
    auto sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    plan_weights.emplace_back(parallel_qv, std::move(weights));
    max_weight = std::max(sum, max_weight);
  }

  /* split candidates of parallelizing qv into buckets by limiting each bucket weight */
  double bucket_weight_limit = max_weight / 3.0;
  double max_bucket_weight = 0, max_single_vertex_weight = 0;
  if (verbosePlannerLog()) {
    LOG(INFO) << "===== Parallelizing plans =====";
    LOG(INFO) << "max weight " << max_weight << ", limit " << bucket_weight_limit;
    LOG(INFO) << "Partitioned plan count before " << partitioned_plans->size();
  }
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> parallelized_partitioned_plans;
  parallelized_partitioned_plans.reserve(partitioned_plans->size());
  QueryVertexID parallelizing_qv = DUMMY_QUERY_VERTEX;
  for (uint32_t i = 0; i < partitioned_plans->size(); ++i) {
    auto& pair = (*partitioned_plans)[i];
    parallelizing_qv = plan_weights[i].first;
    const auto& weights = plan_weights[i].second;

    double bucket_weight = 0;
    for (uint32_t l = 0, r = 0; r < weights.size(); ++r) {
      bucket_weight += weights[r];
      max_single_vertex_weight = std::max(bucket_weight, max_single_vertex_weight);
      if (bucket_weight > bucket_weight_limit) {
        if (l == r) {
          if (bucket_weight > max_bucket_weight) {
            max_bucket_weight = bucket_weight;
          }
          parallelized_partitioned_plans.emplace_back(pair);
          parallelized_partitioned_plans.back().second[parallelizing_qv].addRange(l, r);
          if (verbosePlannerLog()) {
            LOG(INFO) << "Partition " << (parallelized_partitioned_plans.size() - 1) << " plan " << pair.first
                      << " weight " << bucket_weight;
          }
          l = r + 1;
          bucket_weight = 0;
        } else {
          bucket_weight -= weights[r];
          if (bucket_weight > max_bucket_weight) {
            max_bucket_weight = bucket_weight;
          }
          parallelized_partitioned_plans.emplace_back(pair);
          parallelized_partitioned_plans.back().second[parallelizing_qv].addRange(l, r - 1);
          if (verbosePlannerLog()) {
            LOG(INFO) << "Partition " << (parallelized_partitioned_plans.size() - 1) << " plan " << pair.first
                      << " weight " << bucket_weight;
          }
          l = r;
          --r;
          bucket_weight = 0;
        }
      } else if (r == weights.size() - 1) {
        if (bucket_weight > max_bucket_weight) {
          max_bucket_weight = bucket_weight;
        }
        parallelized_partitioned_plans.emplace_back(std::move(pair));
        if (verbosePlannerLog()) {
          LOG(INFO) << "Partition " << (parallelized_partitioned_plans.size() - 1) << " plan " << pair.first
                    << " weight " << bucket_weight;
        }
        if (l > 0) {
          parallelized_partitioned_plans.back().second[parallelizing_qv].addRange(l, r);
        }
      }
    }
  }
  partitioned_plans->swap(parallelized_partitioned_plans);

  if (verbosePlannerLog()) {
    LOG(INFO) << "max_bucket_weight " << max_bucket_weight << " single vertex max weight " << max_single_vertex_weight;
  }
}

// TODO(engineering): classes in plan should not know classes in exec
BacktrackingPlan* Planner::generateExecutionPlan(const CandidateResult* result, bool multithread) {
  backtracking_plan_ = std::make_unique<BacktrackingPlan>();
  if (!toPartitionCandidates()) {  // no partition, generate one logical plan for single-threaded execution
    auto candidate_views = result->getCandidates();
    auto cardinality = getCardinality(candidate_views);
    backtracking_plan_->addPlan(generateLogicalExecutionPlan(cardinality, DUMMY_QUERY_VERTEX, &candidate_views));
    return backtracking_plan_.get();
  }

  auto n_partitions = query_context_->graph_metadata->numPartitions();
  auto partitioned_result = (PartitionedCandidateResult*)result;

  // 1. select partitioning query vertices
  auto partitioning_qv = getPartitioningQueryVertices();
  if (shortPlannerLog()) {
    std::stringstream ss;
    for (auto v : partitioning_qv) {
      ss << ' ' << v;
    }
    LOG(INFO) << partitioning_qv.size() << " partitioning query vertices:" << ss.str();
  }

  // 2. generate scopes of partitioned backtracking search space
  auto partitions = generatePartitionedScopes(n_partitions, partitioning_qv, partitioned_result);

  auto t0 = std::chrono::steady_clock::now();
  // 3. for each scope, generate a logical execution plan
  auto partition_plans = generateLogicalPlans(partitioning_qv, partitions, partitioned_result);

  if (verbosePlannerLog()) {
    for (auto& pair : partition_plans) {
      std::stringstream ss;
      for (auto& scope : pair.second) {
        ss << '\t';
        scope.print(ss);
      }
      LOG(INFO) << "Plan " << pair.first << ss.str();
    }
  }

  auto t1 = std::chrono::steady_clock::now();
  LOG(INFO) << ">>>>>>>>>> Time to generateLogicalPlans " << toSeconds(t0, t1) << "s";
  // 4. parallelize partitioned plans for better concurrency and load balance
  if (multithread) {
    parallelizePartitionedPlans(partitioning_qv, &partition_plans, partitioned_result);
    CHECK_EQ(backtracking_plan_->getPlans().size(), backtracking_plan_->getNumInputOperators());

    auto t2 = std::chrono::steady_clock::now();
    LOG(INFO) << ">>>>>>>>>> Time to parallelizePartitionedPlans " << toSeconds(t1, t2) << "s";
  }

  // TODO(tatiana): check for prunable partitioned plans

  if (shortPlannerLog()) {  // debug log for scopes and plans
    uint32_t plan_idx = 0;
    for (auto plan : backtracking_plan_->getPlans()) {
      LOG(INFO) << "----- Plan " << plan_idx++ << " -----";
      plan->printPhysicalPlan();
    }
    for (auto& pair : partition_plans) {
      std::stringstream ss;
      for (auto& scope : pair.second) {
        ss << '\t';
        scope.print(ss);
      }
      LOG(INFO) << "Plan " << pair.first << ss.str();
    }
  }

  backtracking_plan_->addPartitionedPlans(std::move(partition_plans));
  return backtracking_plan_.get();
}

BacktrackingPlan* Planner::generateExecutionPlan(const std::vector<std::vector<VertexID>>* candidate_cardinality,
                                                 bool multithread) {
  // TODO(tatiana): for explain mode
  return backtracking_plan_.get();
}

ExecutionPlan* Planner::generateLogicalExecutionPlan(const std::vector<VertexID>& candidate_cardinality,
                                                     QueryVertexID require_only,
                                                     const std::vector<CandidateSetView>* candidate_views,
                                                     const std::vector<QueryVertexID>* partitioning_qvs,
                                                     const std::vector<QueryVertexID>* use_order) {
  if (verbosePlannerLog()) {
    std::stringstream ss;
    for (auto c : candidate_cardinality) {
      ss << ' ' << c;
    }
    LOG(INFO) << "Candidate sizes" << ss.str();
  }
  std::vector<double> cardinality{candidate_cardinality.begin(), candidate_cardinality.end()};

  /* if require_only is a valid query vertex, require the candidate set for require_only to be non-empty,
   * but generate an execution plan even if there is an empty candidate set for other query vertices */
  if (require_only != DUMMY_QUERY_VERTEX) {
    CHECK_LT(require_only, cardinality.size());
    if (cardinality[require_only] == 0) return nullptr;
    for (auto& c : cardinality) {
      if (c == 0) {
        c = 2;
      }
    }
  } else {
    for (auto& c : cardinality) {
      if (c == 0) {
        return nullptr;
      }
    }
  }

  /* The data graph representation varies depending on the execution strategy.
   * For query execution with partitioned graphs, use GraphView. For query on normal graphs, use Normal;
   * For query using an auxiliary bipartite-graph-based index, use BipartiteGraphView. */
  GraphType graph_type =
      toPartitionCandidates()
          ? GraphType::GraphView
          : (query_context_->query_config.use_auxiliary_index
                 ? GraphType::BipartiteGraphView
                 : (query_context_->graph_metadata->numPartitions() > 1 ? GraphType::Partitioned : GraphType::Normal));
  planners_.push_back(std::make_unique<NaivePlanner>(&query_context_->query_graph,
                                                     query_context_->query_config.use_two_hop_traversal,
                                                     std::move(cardinality), graph_type));
  auto& planner = planners_.back();

  auto& order = planner->generateOrder(query_context_->data_graph, *query_context_->graph_metadata, candidate_views,
                                       candidate_cardinality, query_context_->query_config.order_strategy, use_order);
  if (order.empty()) {
    return nullptr;
  }

  ExecutionPlan* plan = nullptr;
  switch (query_context_->query_config.compression_strategy) {
  case CompressionStrategy::Static: {
    plan = planner->generatePlan();
    plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()));
    break;
  }
  case CompressionStrategy::Dynamic: {
    plan = planner->generatePlanWithDynamicCover(query_context_->data_graph, candidate_views);
    plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()) &&
                          (plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0));
    break;
  }
  case CompressionStrategy::None: {
    plan = planner->generatePlanWithoutCompression();
    plan->setInputAreKeys(true);
  }
  }

  // one logical input operator for each logical plan
  newInputOperators(plan, partitioning_qvs);
  return plan;
}

}  // namespace circinus

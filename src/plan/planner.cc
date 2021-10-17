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

#include "algorithms/automorphism_check.h"
#include "algorithms/centrality.h"
#include "ops/logical/compressed_input.h"
#include "ops/order_generator.h"
#include "plan/backtracking_plan.h"
#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/query_utils.h"

namespace circinus {

/**
 * For now, in the case of partitioned graph, all graph partitions share the same candidate pruning plan.
 * TODO(tatiana): cache order / bfs-tree for execution plan generation
 */
CandidatePruningPlan* Planner::generateCandidatePruningPlan() {
  auto& q = query_context_->query_graph;
  auto& strategy = query_context_->query_config.candidate_pruning_strategy;
  switch (strategy) {
  case CandidatePruningStrategy::None:
  case CandidatePruningStrategy::Online: {
    if (query_context_->query_config.seed.first != DUMMY_QUERY_VERTEX) {
      candidate_pruning_plan_.setFinished();
    } else {
      OrderGenerator order(&q);
      auto start_qv = order.getOrder(query_context_->query_config.order_strategy, DUMMY_QUERY_VERTEX).front();
      candidate_pruning_plan_.newLDFScan(q, {start_qv});
      candidate_pruning_plan_.setPartitionResult(false);
      // candidate_pruning_plan_.setPartitionResult(toPartitionCandidates());
    }
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
  QueryVertexID seed_qv = query_context_->query_config.seed.first;
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
    // switch (strategy) {
    // case CandidatePruningStrategy::CFL: {
    //   candidate_pruning_plan_.newCFLFilter(&q, metadata, cardinality, seed_qv);
    //   break;
    // }
    // default:
    //   candidate_pruning_plan_.setFinished();
    // }
    // return &candidate_pruning_plan_;
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
      candidate_pruning_plan_.newCFLFilter(&q, metadata, cardinality, seed_qv);
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

std::pair<int, std::vector<std::pair<uint32_t, std::vector<CandidateScope>>>> Planner::generateLogicalPlans(
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
        if (shortPlannerLog()) {
          LOG(INFO) << "Generate Plan " << backtracking_plan_->getPlans().size() << " for intra-partition "
                    << graph_partition;
        }
        std::vector<CandidateScope> plan_scopes(n_qvs, scopes[primary_pqv]);
        auto& stats = candidate_cardinality[graph_partition];
        auto candidate_views = result->getCandidatesByScopes(plan_scopes);
        auto plan = generateLogicalExecutionPlan(stats, primary_pqv, &candidate_views, &partitioning_qvs);
        if (plan == nullptr) {
          if (planners_.size() == backtracking_plan_->getPlans().size() + 1) {  // need a global plan
            planners_.pop_back();
            if (global_plan_idx == -1) {
              auto candidate_views = result->getCandidates();
              // for global plan use candidate cardinality product for lightweight cover selection
              plan = generateLogicalExecutionPlan(getCardinality(candidate_views), DUMMY_QUERY_VERTEX, &candidate_views,
                                                  &partitioning_qvs, nullptr, false);
              CHECK(plan != nullptr);  // the entire problem does not have a match if nullptr
              plan->setStepCosts(std::vector<double>(plan->getStepCosts().size()));
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
    if (global_plan_idx != -1) {
      LOG(INFO) << "global plan idx " << global_plan_idx;
    }
    return {global_plan_idx, partition_plans};
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
  return {-1, partition_plans};
}

void Planner::parallelizePartitionedPlans(
    const std::vector<QueryVertexID>& partitioning_qvs,
    std::vector<std::pair<uint32_t, std::vector<CandidateScope>>>* partitioned_plans,
    PartitionedCandidateResult* result) {
  parallelizePartitionedPlans(partitioning_qvs, partitioned_plans, result, 5e5, 5e5);
}

void Planner::parallelizePartitionedPlans(
    const std::vector<QueryVertexID>& partitioning_qvs,
    std::vector<std::pair<uint32_t, std::vector<CandidateScope>>>* partitioned_plans,
    PartitionedCandidateResult* result, double bucket_weight_limit, double parallelization_threshold) {
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> parallelized_partitioned_plans;
  std::vector<bool> segment_mask;
  parallelized_partitioned_plans.reserve(partitioned_plans->size());
  segment_mask.reserve(partitioned_plans->size());
  for (auto& partition : *partitioned_plans) {
    auto& step_costs = backtracking_plan_->getPlan(partition.first)->getStepCosts();
    auto sum_costs = std::accumulate(step_costs.begin(), step_costs.end(), 0.) + 1.;
    if (sum_costs <= parallelization_threshold) {
      parallelized_partitioned_plans.push_back(std::move(partition));
      segment_mask.push_back(false);
      continue;
    }

    auto candidate_views = result->getCandidatesByScopes(partition.second);
    auto& plan = backtracking_plan_->getPlan(partition.first);
    auto parallel_qv = plan->getMatchingOrder().front();
    if (!plan->isInCover(parallel_qv)) {
      parallel_qv = plan->getMatchingOrder()[1];
      CHECK(plan->isInCover(parallel_qv));
    }

    auto weights = planners_[partition.first]->getParallelizingQueryVertexWeights(
        parallel_qv, query_context_->data_graph, candidate_views, plan->getQueryCoverBits());
    auto sum_weight = std::accumulate(weights.begin(), weights.end(), 0);

    // bucket_weight_limit = sum_weight * parallelization_threshold / sum_costs;
    auto partition_bucket_weight_limit = std::max(bucket_weight_limit, sum_weight / 100.);
    if (sum_weight < partition_bucket_weight_limit) {
      parallelized_partitioned_plans.push_back(std::move(partition));
      segment_mask.push_back(true);
      continue;
    }
    double bucket_weight = 0;
    for (uint32_t l = 0, r = 0; r < weights.size(); ++r) {  // l and r both inclusive
      bucket_weight += weights[r];
      if (bucket_weight > partition_bucket_weight_limit) {
        if (l != r) {
          bucket_weight -= weights[r];
          --r;
        }
        parallelized_partitioned_plans.emplace_back(partition);
        parallelized_partitioned_plans.back().second[parallel_qv].addRange(l, r);
        if (verbosePlannerLog()) {
          LOG(INFO) << "Partition " << (parallelized_partitioned_plans.size() - 1) << " plan " << partition.first
                    << " weight " << bucket_weight;
        }
        l = r + 1;
        bucket_weight = 0;
        segment_mask.push_back(true);
      } else if (r == weights.size() - 1) {
        parallelized_partitioned_plans.emplace_back(std::move(partition));
        if (verbosePlannerLog()) {
          LOG(INFO) << "Partition " << (parallelized_partitioned_plans.size() - 1) << " plan " << partition.first
                    << " weight " << bucket_weight;
        }
        CHECK_GT(l, 0);
        parallelized_partitioned_plans.back().second[parallel_qv].addRange(l, r);
        segment_mask.push_back(true);
      }
    }
  }

  DCHECK_EQ(segment_mask.size(), parallelized_partitioned_plans.size());
  partitioned_plans->swap(parallelized_partitioned_plans);
  backtracking_plan_->setPartitionedPlanSegmentMask(std::move(segment_mask));
}

std::vector<std::pair<uint32_t, std::vector<QueryVertexID>>> Planner::generateParallelQueryVertex(
    const std::vector<std::pair<uint32_t, std::vector<CandidateScope>>>& partition_plans) {
  std::vector<std::pair<uint32_t, std::vector<QueryVertexID>>> parallel_opids(partition_plans.size());
  for (uint32_t pidx = 0; pidx < partition_plans.size(); ++pidx) {
    auto& partition = partition_plans[pidx];
    QueryVertexID parallel_qv = planners_[partition.first]->selectParallelizingQueryVertex(
        backtracking_plan_->getPlan(partition.first)->getQueryCoverBits(), std::vector<QueryVertexID>{});
    if (parallel_qv == DUMMY_QUERY_VERTEX) {
      parallel_opids[pidx] = std::make_pair(partition.first, std::vector<QueryVertexID>());
      continue;
    }
    const std::vector<QueryVertexID>& matching_order = planners_[partition.first]->getMatchingOrder();
    // find matching order of parallel_qv
    uint32_t parallel_qv_order =
        std::distance(matching_order.begin(), std::find(matching_order.begin(), matching_order.end(), parallel_qv));
    // find the operator index corresponding to extpanding parallel_qv
    std::vector<Operator*>& ops = backtracking_plan_->getOperators(partition.first);
    uint32_t prefix_vertex_num = 0;
    for (uint32_t i = 0; i < ops.size() - 1; ++i) {
      if (dynamic_cast<TraverseOperator*>(ops[i])->extend_vertex()) {
        ++prefix_vertex_num;
      }
      if (prefix_vertex_num == parallel_qv_order) {
        std::vector<QueryVertexID> opid{i + 1};
        parallel_opids[pidx] = std::make_pair(partition.first, std::move(opid));
        break;
      }
    }
  }
  return std::move(parallel_opids);
}

BacktrackingPlan* Planner::generateExecutionPlan(std::pair<QueryVertexID, VertexID> seed, bool multithread) {
  backtracking_plan_ = std::make_unique<BacktrackingPlan>();
  const QueryGraph& query_graph = query_context_->query_graph;
  const ReorderedPartitionedGraph* data_graph = dynamic_cast<ReorderedPartitionedGraph*>(query_context_->data_graph);
  if (query_graph.getVertexLabel(seed.first) != ALL_LABEL &&
      query_graph.getVertexLabel(seed.first) != data_graph->getVertexLabel(seed.second)) {
    return nullptr;
  }

  // FIXME(tatiana): seed qv should be excluded from symmetry analysis breakSymmetry();
  planners_.push_back(std::make_unique<NaivePlanner>(&query_context_->query_graph, getGraphType()));

  auto& planner = planners_.back();

  auto& order = planner->generateOrder(seed.first, query_context_->query_config.order_strategy, qv_partial_order_);

  if (order.empty()) {
    return nullptr;
  }

  if (shortPlannerLog()) {
    LOG(INFO) << "Order:" << toString(order);
  }

  ExecutionPlan* plan = planner->generatePlanWithDynamicCover(query_context_->data_graph, nullptr, qv_partial_order_);
  plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()));

  backtracking_plan_->addPlan(plan);
  return backtracking_plan_.get();
}

void Planner::breakSymmetry() {
  DCHECK(backtracking_plan_ != nullptr);
  LOG(INFO) << "Symmetry breaking enabled? " << FLAGS_break_symmetry;
  if (FLAGS_break_symmetry) {
    AutomorphismCheck ac(query_context_->query_graph);
    auto conds = ac.getPartialOrder();
    for (auto& pair : conds) {
      qv_partial_order_[pair.first].second.push_back(pair.second);
      qv_partial_order_[pair.second].first.push_back(pair.first);
    }
    if (shortPlannerLog()) {
      std::stringstream ss;
      for (auto& pair : conds) {
        ss << ' ' << pair.first << '<' << pair.second;
      }
      LOG(INFO) << "Partial orders:" << ss.str();
    }
    // TODO(tatiana): populate the conditions by transtivity for probably earlier pruning?
    backtracking_plan_->setQueryPartialOrder(qv_partial_order_);
  }
}

// TODO(engineering): classes in plan should not know classes in exec
BacktrackingPlan* Planner::generateExecutionPlan(const CandidateResult* result, bool multithread) {
  backtracking_plan_ = std::make_unique<BacktrackingPlan>();
  breakSymmetry();
  if (query_context_->query_config.candidate_pruning_strategy ==
      CandidatePruningStrategy::Online) {  // fast plan: only one plan, no partitioning applied
    // generate order and compression plan
    planners_.push_back(std::make_unique<NaivePlanner>(&query_context_->query_graph, getGraphType()));
    auto& planner = planners_.back();
    auto& order = planner->generateOrder(query_context_->query_config.seed.first,
                                         query_context_->query_config.order_strategy, qv_partial_order_);

    if (shortPlannerLog()) {
      LOG(INFO) << "Order:" << toString(order);
    }
    ExecutionPlan* plan = nullptr;
    switch (query_context_->query_config.compression_strategy) {
    case CompressionStrategy::Static: {
      plan = planner->generatePlan(qv_partial_order_);
      plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()));
      break;
    }
    case CompressionStrategy::Dynamic: {
      plan = planner->generatePlanWithDynamicCover(query_context_->data_graph, nullptr, qv_partial_order_);
      plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()) &&
                            (plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0));
      break;
    }
    case CompressionStrategy::None: {
      plan = planner->generatePlanWithoutCompression();
      plan->setInputAreKeys(true);
    }
    }

    // parallelize plan
    if (multithread) {
      VertexID start_size = 0;
      if (result == nullptr) {
        start_size = query_context_->graph_metadata->getNumVertices();
      } else {
        auto candidates = result->getCandidates();
        CHECK_EQ(candidates.size(), 1);
        start_size = candidates.front().size();
      }
      const uint32_t n_chunks = 100;
      std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> partitioned_plans;
      if (start_size < n_chunks) {
        partitioned_plans.resize(start_size);
        for (uint32_t i = 0; i < start_size; ++i) {
          partitioned_plans[i].first = 0;
          CandidateScope scope;
          partitioned_plans[i].second.resize(1);
          partitioned_plans[i].second.back().addRange(i, i);  // both inclusive
        }
      } else {
        partitioned_plans.resize(n_chunks);
        uint32_t chunk_size = start_size / n_chunks;
        auto extra = start_size % n_chunks;
        for (uint32_t i = 0; i < n_chunks; ++i) {
          partitioned_plans[i].first = 0;
          CandidateScope scope;
          partitioned_plans[i].second.resize(1);
          if (i < extra) {
            partitioned_plans[i].second.back().addRange(i * (chunk_size + 1),
                                                        (i + 1) * (chunk_size + 1) - 1);  // both inclusive
          } else {
            partitioned_plans[i].second.back().addRange(i * chunk_size + extra,
                                                        (i + 1) * chunk_size + extra - 1);  // both inclusive
          }
        }
      }
      backtracking_plan_->addPartitionedPlans(std::move(partitioned_plans));
    } else {
      std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> partitioned_plans(1);
      partitioned_plans.front().first = 0;
      partitioned_plans.front().second.resize(1);  // full scope
      backtracking_plan_->addPartitionedPlans(std::move(partitioned_plans));
    }

    backtracking_plan_->addInputOperator(std::make_unique<LogicalCompressedInputOperator>(0, plan->inputAreKeys()));
    backtracking_plan_->setPartitionedPlanSegmentMask(
        std::vector<bool>(backtracking_plan_->getNumPartitionedPlans(), true));
    backtracking_plan_->addPlan(plan);
    return backtracking_plan_.get();
  }

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
  auto[global_plan_idx, partition_plans] = generateLogicalPlans(partitioning_qv, partitions, partitioned_result);

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

  CHECK_GT(partition_plans.size(), 0);
  backtracking_plan_->addPartitionedPlans(std::move(partition_plans));

  return backtracking_plan_.get();
}

BacktrackingPlan* Planner::generateExecutionPlan(const std::vector<std::vector<VertexID>>* candidate_cardinality,
                                                 bool multithread) {
  // TODO(engineering): for explain mode
  return backtracking_plan_.get();
}

ExecutionPlan* Planner::generateLogicalExecutionPlan(const std::vector<VertexID>& candidate_cardinality,
                                                     QueryVertexID require_only,
                                                     const std::vector<CandidateSetView>* candidate_views,
                                                     const std::vector<QueryVertexID>* partitioning_qvs,
                                                     const std::vector<QueryVertexID>* use_order, bool use_cover_path) {
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

  GraphType graph_type = toPartitionCandidates() ? GraphType::GraphView : getGraphType();
  planners_.push_back(std::make_unique<NaivePlanner>(&query_context_->query_graph,
                                                     query_context_->query_config.use_two_hop_traversal,
                                                     std::move(cardinality), graph_type));
  auto& planner = planners_.back();

  auto& order = planner->generateOrder(query_context_->data_graph, *query_context_->graph_metadata, candidate_views,
                                       candidate_cardinality, query_context_->query_config.order_strategy,
                                       query_context_->query_config.seed.first, use_order);

  if (order.empty()) {
    return nullptr;
  }

  if (shortPlannerLog()) {
    LOG(INFO) << "Order:" << toString(order);
  }

  ExecutionPlan* plan = nullptr;
  switch (query_context_->query_config.compression_strategy) {
  case CompressionStrategy::Static: {
    plan = planner->generatePlan();
    plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()));
    break;
  }
  case CompressionStrategy::Dynamic: {
    plan =
        planner->generatePlanWithDynamicCover(query_context_->data_graph, use_cover_path ? candidate_views : nullptr);
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

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
#include <utility>
#include <vector>

#include "algorithms/centrality.h"
#include "ops/logical/compressed_input.h"
#include "ops/order_generator.h"
#include "plan/backtracking_plan.h"
#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
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

// FIXME(tatiana): classes in plan should not know classes in exec
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
  std::stringstream ss;
  for (auto v : order_by_closeness_centrality) {
    ss << ' ' << v << ": " << centrality[v] << ",  ";
  }
  LOG(INFO) << "Query vertex centrality: [" << ss.str() << ']';

  order_by_closeness_centrality.resize(1);

  return order_by_closeness_centrality;
}

void Planner::exhaustivePartitionPlan(std::vector<std::vector<CandidateScope>>& partitioned_plans, uint32_t level,
                                      uint32_t partition_num, const std::vector<QueryVertexID>& partitioning_qv,
                                      std::vector<CandidateScope> partitioned_plan) {
  if (level == partitioning_qv.size()) {
    std::stringstream ss;
    for (auto x : partitioned_plan) {
      x.print(ss);
    }
    LOG(INFO) << ss.str();
    partitioned_plans.emplace_back(std::move(partitioned_plan));
    return;
  }

  for (uint32_t i = 0; i < partition_num; ++i) {
    std::vector<CandidateScope> new_partitioned_plan(partitioned_plan.begin(), partitioned_plan.end());
    new_partitioned_plan[partitioning_qv[level]].usePartition(i);
    exhaustivePartitionPlan(partitioned_plans, level + 1, partition_num, partitioning_qv, new_partitioned_plan);
  }
}

std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> Planner::generatePartitionedPlans(
    const std::vector<QueryVertexID>& partitioning_qv) {
  // CHECK_EQ(partitioning_qv.size(), 1);  // consider the locality indicator only, no. tasks equal to no. partitions
  DLOG(INFO) << "partitioning_qv = " << partitioning_qv.front();
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> ret;
  ret.reserve(backtracking_plan_->getPlans().size() * ((1 << partitioning_qv.size()) - 1));

  auto n_qvs = query_context_->query_graph.getNumVertices();

  for (uint32_t i = 0; i < backtracking_plan_->getPlans().size(); ++i) {
    ret.emplace_back(i, std::vector<CandidateScope>(n_qvs));
    // for (uint32_t j = 0; j < n_qvs; ++j)
    ret.back().second[partitioning_qv.front()].usePartition(backtracking_plan_->getPlan(i)->getPartitionId());
    // for (uint64_t j = 1; j < all; ++j) {
    //   uint64_t tmp = j;
    //   ret.emplace_back(i, std::vector<CandidateScope>(n_qvs));
    //   uint32_t idx = 0;
    //   while (tmp) {
    //     if (tmp % 3 == 0) {
    //       ret.back().second[partitioning_qv[idx]].beforePartition(backtracking_plan_->getPlan(i)->getPartitionId());
    //     } else if (tmp % 3 == 1) {
    //       ret.back().second[partitioning_qv[idx]].usePartition(backtracking_plan_->getPlan(i)->getPartitionId());
    //     } else if (tmp % 3 == 2) {
    //       ret.back().second[partitioning_qv[idx]].afterPartition(backtracking_plan_->getPlan(i)->getPartitionId());
    //     }
    //     idx++;
    //   }
    // }
  }
  LOG(INFO) << "Partition qeury vertex size: " << partitioning_qv.size() << ",  partition plan size:" << ret.size();
  return ret;
}

void Planner::newInputOperators(const QueryGraph& q, const std::vector<QueryVertexID>& partitioning_qvs) {
  for (auto logical_plan : backtracking_plan_->getPlans()) {
    // TODO(byli): use neighborhood filter corresponding to the candidate pruning strategy?
    backtracking_plan_->addInputOperator(std::make_unique<PartitionedLogicalCompressedInputOperator>(
        &q, logical_plan->inputAreKeys(), logical_plan->getMatchingOrder(), partitioning_qvs));
  }
}

void Planner::newInputOperators() {
  auto& logical_plan = backtracking_plan_->getPlans().front();
  backtracking_plan_->addInputOperator(std::make_unique<LogicalCompressedInputOperator>(
      logical_plan->getMatchingOrder().front(), logical_plan->inputAreKeys()));
}

// FIXME(tatiana): classes in plan should not know classes in exec
BacktrackingPlan* Planner::generateExecutionPlan(const CandidateResult* result, bool multithread) {
  auto candidate_cardinality = result->getCandidateCardinality();
  if (!toPartitionCandidates()) {  // no partition, generate one execution plan
    auto candidate_views = result->getCandidates();
    backtracking_plan_ = std::make_unique<BacktrackingPlan>();
    if (candidate_cardinality.size() > 1) {
      std::vector<VertexID> sum(candidate_cardinality.front().size(), 0);
      for (auto& part : candidate_cardinality) {
        for (uint32_t i = 0; i < part.size(); ++i) {
          sum[i] += part[i];
        }
      }
      auto order_generator = OrderGenerator(query_context_->data_graph, *query_context_->graph_metadata,
                                            &query_context_->query_graph, candidate_views, sum);
      auto use_order = order_generator.getOrder(query_context_->query_config.order_strategy);
      backtracking_plan_->addPlan(generateExecutionPlan(sum, use_order, multithread, &candidate_views));
    } else {
      auto order_generator =
          OrderGenerator(query_context_->data_graph, *query_context_->graph_metadata, &query_context_->query_graph,
                         candidate_views, candidate_cardinality.front());
      auto use_order = order_generator.getOrder(query_context_->query_config.order_strategy);
      backtracking_plan_->addPlan(
          generateExecutionPlan(candidate_cardinality.front(), use_order, multithread, &candidate_views));
    }
    // one logical input operator
    newInputOperators();
    return backtracking_plan_.get();
  }

  backtracking_plan_ = std::make_unique<BacktrackingPlan>();
  // generate a logical plan for each partition
  uint32_t n_qvs = query_context_->query_graph.getNumVertices();
  result = (PartitionedCandidateResult*)result;
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> partition_plans;
  if (false) {
    auto partitioning_qv = getPartitioningQueryVerticesByClosenessCentrality();
    std::vector<std::vector<CandidateScope>> partitioned_plans;
    exhaustivePartitionPlan(partitioned_plans, 0, candidate_cardinality.size(), partitioning_qv,
                            std::vector<CandidateScope>(n_qvs));
    partition_plans.reserve(partitioned_plans.size());
    uint32_t partition_plan_num = 0;
    for (auto& scopes : partitioned_plans) {
      LOG(INFO) << "Generate Plan " << partition_plan_num;
      std::vector<uint64_t> stats(n_qvs);
      auto candidate_views = ((PartitionedCandidateResult*)result)->getCandidatesByScopes(scopes);
      for (uint32_t i = 0; i < n_qvs; ++i) {
        stats[i] = candidate_views[i].size();
      }
      auto order_generator = OrderGenerator(query_context_->data_graph, *query_context_->graph_metadata,
                                            &query_context_->query_graph, candidate_views, stats);
      auto use_order = order_generator.getOrder(query_context_->query_config.order_strategy);
      auto plan = generateExecutionPlan(stats, use_order, multithread, &candidate_views, partitioning_qv.front());
      if (plan == nullptr) {
        continue;
      }
      backtracking_plan_->addPlan(plan);
      partition_plans.emplace_back(partition_plan_num++, scopes);
    }

    backtracking_plan_->addPartitionedPlans(std::move(partition_plans));
    // one logical input operator for each logical plan
    LOG(INFO) << backtracking_plan_->getPartitionedPlan(0).second.size();
    newInputOperators(query_context_->query_graph, partitioning_qv);
  } else {
    LOG(INFO) << "-----------------";
    auto partitioning_qv = getPartitioningQueryVerticesByClosenessCentrality();
    double max_weight = 0;
    for (uint32_t i = 0; i < candidate_cardinality.size(); ++i) {
      LOG(INFO) << "Generate Plan " << i;
      // TODO(tatiana): apply the state-of-the-art matching order strategies, compute the orders only
      std::vector<CandidateScope> scopes(n_qvs);
      auto& stats = candidate_cardinality[i];
      for (auto& scope : scopes) {
        scope.usePartition(i);
      }
      auto candidate_views = ((PartitionedCandidateResult*)result)->getCandidatesByScopes(scopes);

      auto order_generator = OrderGenerator(query_context_->data_graph, query_context_->graph_metadata->getPartition(i),
                                            &query_context_->query_graph, candidate_views, stats);
      auto use_order = order_generator.getOrder(query_context_->query_config.order_strategy);
      // use_order = {0, 2, 1, 6, 3, 7, 4, 5};
      // use_order = {3, 0, 1, 2, 7, 4, 6, 5};
      // use_order = {1, 2, 0, 3, 5, 4, 6, 7};
      auto plan = generateExecutionPlan(stats, use_order, multithread, &candidate_views, partitioning_qv.front());
      backtracking_plan_->addPlan(plan);
      if (plan->getPartitionQueryVertexWeightsSum() > max_weight) {
        max_weight = plan->getPartitionQueryVertexWeightsSum();
      }
    }

    double bucket_weight_limit = max_weight / 3.0;
    // double bucket_weight_limit = 0;
    double x = 0, y = 0;
    LOG(INFO) << "max weight " << max_weight << ", limit " << bucket_weight_limit;
    for (uint32_t i = 0; i < backtracking_plan_->getPlans().size(); ++i) {
      auto plan = backtracking_plan_->getPlan(i);
      const auto& weights = plan->getPartitionQueryVertexWeights();
      // partition_plans.emplace_back(i, std::vector<CandidateScope>(n_qvs));
      // partition_plans.back().second[partitioning_qv.front()].addRange(i, 1, 1);
      // partition_plans.back().second[partitioning_qv.front()].usePartition(i);
      LOG(INFO) << "----------------------" << 0 << " " << weights.size() - 1;
      double bucket_weight = 0;
      for (uint32_t l = 0, r = 0; r < weights.size(); ++r) {
        bucket_weight += weights[r];
        if (bucket_weight > bucket_weight_limit) {
          if (l == r) {
            if (bucket_weight > x) {
              x = bucket_weight;
              y = x;
            }
            partition_plans.emplace_back(i, std::vector<CandidateScope>(n_qvs));
            partition_plans.back().second[partitioning_qv.front()].addRange(i, l, r);
            LOG(INFO) << l << " " << r;
            l = r + 1;
            bucket_weight = 0;
          } else {
            if (bucket_weight - weights[r] > x) {
              x = bucket_weight - weights[r];
            }
            partition_plans.emplace_back(i, std::vector<CandidateScope>(n_qvs));
            partition_plans.back().second[partitioning_qv.front()].addRange(i, l, r - 1);
            LOG(INFO) << l << " " << r - 1;
            l = r;
            --r;
            bucket_weight = 0;
          }
        } else if (r == weights.size() - 1) {
          if (bucket_weight > x) {
            x = bucket_weight;
          }
          partition_plans.emplace_back(i, std::vector<CandidateScope>(n_qvs));
          partition_plans.back().second[partitioning_qv.front()].addRange(i, l, r);
          LOG(INFO) << l << " " << r;
        }
      }
    }

    LOG(INFO) << "max_bucket_weight " << x << " single vertex max weight " << y;
    backtracking_plan_->addPartitionedPlans(std::move(partition_plans));
    // one logical input operator for each logical plan
    newInputOperators(query_context_->query_graph, partitioning_qv);
  }
  LOG(INFO) << "Generated plan";
  return backtracking_plan_.get();
}

BacktrackingPlan* Planner::generateExecutionPlan(const std::vector<std::vector<VertexID>>* candidate_cardinality,
                                                 bool multithread) {
  auto use_order = getOrder(query_context_->query_config.matching_order, query_context_->query_graph.getNumVertices());
  if (!toPartitionCandidates()) {  // no partition, generate one execution plan
    backtracking_plan_ = std::make_unique<BacktrackingPlan>();
    if (candidate_cardinality->size() > 1) {
      std::vector<VertexID> sum(candidate_cardinality->front().size(), 0);
      for (auto& part : *candidate_cardinality) {
        for (uint32_t i = 0; i < part.size(); ++i) {
          sum[i] += part[i];
        }
      }
      backtracking_plan_->addPlan(generateExecutionPlan(sum, use_order, multithread));
    } else {
      auto car = candidate_cardinality->front();
      backtracking_plan_->addPlan(generateExecutionPlan(car, use_order, multithread));
    }
    // one logical input operator
    newInputOperators();
    return backtracking_plan_.get();
  }

  backtracking_plan_ = std::make_unique<BacktrackingPlan>();
  // generate a logical plan for each partition
  for (uint32_t i = 0; i < candidate_cardinality->size(); ++i) {
    auto stats = (*candidate_cardinality)[i];
    auto plan = generateExecutionPlan(stats, use_order, multithread);
    if (plan == nullptr) {
      continue;
    }
    backtracking_plan_->addPlan(plan);
  }
  auto partitioning_qv = getPartitioningQueryVertices();
  backtracking_plan_->addPartitionedPlans(generatePartitionedPlans(partitioning_qv));

  // one logical input operator for each logical plan
  newInputOperators(query_context_->query_graph, partitioning_qv);
  LOG(INFO) << "Generated plan";
  return backtracking_plan_.get();
}

ExecutionPlan* Planner::generateExecutionPlan(std::vector<VertexID>& candidate_cardinality,
                                              const std::vector<QueryVertexID>& use_order, bool multithread,
                                              const std::vector<CandidateSetView>* candidate_views,
                                              QueryVertexID partitioning_qv) {
  // TODO(tatiana): handle the case when the cardinality of a candidate set is 0
  std::stringstream ss;
  for (auto c : candidate_cardinality) {
    ss << c << " ";
  }
  LOG(INFO) << ss.str();
  for (auto& c : candidate_cardinality) {
    if (c == 0) {
      c = 1;
    }
  }
  std::vector<double> cardinality{candidate_cardinality.begin(), candidate_cardinality.end()};

  /* The data graph representation varies depending on the execution strategy.
   * For query execution with partitioned graphs, use GraphView. For query on normal graphs, use Normal;
   * For query using an auxiliary bipartite-graph-based index, use BipartiteGraphView. */
  GraphType graph_type =
      toPartitionCandidates()
          ? GraphType::GraphView
          : (query_context_->query_config.use_auxiliary_index
                 ? GraphType::BipartiteGraphView
                 : (query_context_->graph_metadata->numPartitions() > 1 ? GraphType::Partitioned : GraphType::Normal));
  planners_.push_back(std::make_unique<NaivePlanner>(
      &query_context_->query_graph, query_context_->query_config.use_two_hop_traversal, &cardinality, graph_type));
  auto& planner = planners_.back();

  ExecutionPlan* plan = nullptr;
  // TODO(tatiana): consider make the first vertex in the matching order to be a key if multithread?
  switch (query_context_->query_config.compression_strategy) {
  case CompressionStrategy::Static: {
    plan = planner->generatePlan(use_order);
    plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()));
    break;
  }
  case CompressionStrategy::Dynamic: {
    planner->generateOrder(use_order);
    std::vector<std::vector<double>> cardinality_per_level(cardinality.size(), cardinality);
    planner->generateCoverNode(cardinality_per_level);
    plan = planner->generatePlanWithDynamicCover(query_context_->data_graph, candidate_views, partitioning_qv);
    plan->setInputAreKeys(plan->isInCover(plan->getRootQueryVertexID()) &&
                          (plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0));
    break;
  }
  case CompressionStrategy::None: {
    plan = planner->generatePlanWithoutCompression(use_order);
    plan->setInputAreKeys(true);
  }
  }

  // plan->printPhysicalPlan();
  return plan;
}

}  // namespace circinus

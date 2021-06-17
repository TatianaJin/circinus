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

#include <memory>
#include <string>
#include <vector>

#include "ops/logical/compressed_input.h"
#include "ops/order/cfl_order.h"
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
    candidate_pruning_plan_.setPartitionResult(query_context_->graph_metadata->numPartitions() > 1);
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::NLF: {
    LOG(INFO) << "Prune candidate by NLF";
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
    candidate_pruning_plan_.setPartitionResult(query_context_->graph_metadata->numPartitions() > 1);
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

CandidatePruningPlan* Planner::updateCandidatePruningPlan(const std::vector<std::vector<VertexID>>* part_cardinality) {
  auto& q = query_context_->query_graph;
  auto& strategy = query_context_->query_config.candidate_pruning_strategy;
  auto n_qvs = (*part_cardinality)[0].size();
  std::vector<VertexID> cardinality(n_qvs, 0);
  for (uint32_t i = 0; i < part_cardinality->size(); ++i) {
    for (uint32_t j = 0; j < n_qvs; ++j) {
      cardinality[j] += (*part_cardinality)[i][j];
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
      candidate_pruning_plan_.newDPISOFilter(&q, metadata, cardinality);
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
  candidate_pruning_plan_.setPartitionResult(query_context_->graph_metadata->numPartitions() > 1);
  for (uint32_t i = 0; i < cardinality.size(); ++i) {
    DLOG(INFO) << " |C(v" << i << ")|: " << cardinality[i];
  }
  candidate_pruning_plan_.setFinished();
  return &candidate_pruning_plan_;
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
              return qv_frequency_in_cover[a] >= qv_frequency_in_cover[b];
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

std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> Planner::generatePartitionedPlans(
    const std::vector<QueryVertexID>& partitioning_qv) {
  CHECK_EQ(partitioning_qv.size(), 1);  // consider the locality indicator only, no. tasks equal to no. partitions
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> ret;
  ret.reserve(backtracking_plan_->getPlans().size());

  auto n_qvs = query_context_->query_graph.getNumVertices();
  for (uint32_t i = 0; i < backtracking_plan_->getPlans().size(); ++i) {
    ret.emplace_back(i, std::vector<CandidateScope>(n_qvs));
    ret.back().second[partitioning_qv[0]].usePartition(i);
  }
  return ret;
}

void Planner::newInputOperators(const QueryGraph& q, const GraphMetadata& metadata,
                                const std::vector<std::vector<VertexID>>* candidate_cardinality,
                                const std::vector<QueryVertexID>& partitioning_qvs) {
  for (auto logical_plan : backtracking_plan_->getPlans()) {
    // TODO(byli): use neighborhood filter corresponding to the candidate pruning strategy?
    backtracking_plan_->addInputOperator(std::make_unique<LogicalCompressedInputOperator>(
        logical_plan->inputAreKeys(), logical_plan->getMatchingOrder(), partitioning_qvs));
  }
}

BacktrackingPlan* Planner::generateExecutionPlan(const std::vector<std::vector<VertexID>>* candidate_cardinality,
                                                 bool multithread) {
  if (candidate_cardinality->size() == 1) {  // no partition, generate one execution plan
    backtracking_plan_ = std::make_unique<BacktrackingPlan>();
    backtracking_plan_->addPlan(generateExecutionPlan(candidate_cardinality->front(), multithread));
    return backtracking_plan_.get();
  }

  backtracking_plan_ = std::make_unique<BacktrackingPlan>();
  // generate a logical plan for each partition
  for (auto& stats : *candidate_cardinality) {
    // TODO(tatiana): apply the state-of-the-art matching order strategies, compute the orders only
    backtracking_plan_->addPlan(generateExecutionPlan(stats, multithread));
  }
  auto partitioning_qv = getPartitioningQueryVertices();
  backtracking_plan_->addPartitionedPlans(generatePartitionedPlans(partitioning_qv));

  // one logical input operator for each logical plan
  newInputOperators(query_context_->query_graph, *query_context_->graph_metadata, candidate_cardinality,
                    partitioning_qv);
  LOG(INFO) << "Generated plan";
  return backtracking_plan_.get();
}

ExecutionPlan* Planner::generateExecutionPlan(const std::vector<VertexID>& candidate_cardinality, bool multithread) {
  // TODO(tatiana): handle the case when the cardinality of a candidate set is 0
  for (auto c : candidate_cardinality) {
    CHECK_NE(c, 0);
  }
  std::vector<double> cardinality{candidate_cardinality.begin(), candidate_cardinality.end()};

  /* The data graph representation varies depending on the execution strategy.
   * For query execution with partitioned graphs, use GraphView. For query on normal graphs, use Normal;
   * For query using an auxiliary bipartite-graph-based index, use BipartiteGraphView. */
  GraphType graph_type =
      query_context_->graph_metadata->numPartitions() > 1
          ? GraphType::GraphView
          : (query_context_->query_config.use_auxiliary_index ? GraphType::BipartiteGraphView : GraphType::Normal);
  planners_.push_back(std::make_unique<NaivePlanner>(&query_context_->query_graph, &cardinality, graph_type));
  auto& planner = planners_.back();

  // now order is directly configured
  auto use_order = getOrder(query_context_->query_config.matching_order, query_context_->query_graph.getNumVertices());

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
    plan = planner->generatePlanWithDynamicCover();
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

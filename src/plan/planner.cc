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

#include "ops/order/cfl_order.h"
#include "plan/backtracking_plan.h"
#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/query_utils.h"

namespace circinus {

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
    return &candidate_pruning_plan_;
  }
  case CandidatePruningStrategy::NLF: {
    LOG(INFO) << "Prune candidate by NLF";
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
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

CandidatePruningPlan* Planner::updateCandidatePruningPlan(const std::vector<VertexID>* cardinality) {
  auto& q = query_context_->query_graph;
  auto& strategy = query_context_->query_config.candidate_pruning_strategy;
  if (strategy == CandidatePruningStrategy::LDF || strategy == CandidatePruningStrategy::NLF) {
    for (uint32_t i = 0; i < cardinality->size(); ++i) {
      DLOG(INFO) << "|C(v" << i << ")|: " << (*cardinality)[i];
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
    // TODO(boyang): generate logical neighbor filter for CFL, DAF, TSO and GQL
    auto& metadata = *query_context_->graph_metadata;
    switch (strategy) {
    case CandidatePruningStrategy::DAF: {
      candidate_pruning_plan_.newDPISOFilter(&q, metadata, *cardinality);
      break;
    }
    case CandidatePruningStrategy::CFL: {
      candidate_pruning_plan_.newCFLFilter(&q, metadata, *cardinality);
      break;
    }
    case CandidatePruningStrategy::GQL: {
      candidate_pruning_plan_.newGQLFilter(&q);
      break;
    }
    case CandidatePruningStrategy::TSO: {
      candidate_pruning_plan_.newTSOFilter(&q, metadata, *cardinality);
      break;
    }
    default:
      candidate_pruning_plan_.setFinished();
    }
    return &candidate_pruning_plan_;
  }
  candidate_pruning_plan_.setFinished();
  return &candidate_pruning_plan_;
}

BacktrackingPlan* Planner::generateExecutionPlan(const std::vector<VertexID>* candidate_cardinality, bool multithread) {
  std::vector<double> cardinality{candidate_cardinality->begin(), candidate_cardinality->end()};

  /* The data graph representation varies depending on the execution strategy.
   * For query execution with partitioned graphs, use GraphView. For query on normal graphs, use Normal;
   * For query using an auxiliary bipartite-graph-based index, use BipartiteGraphView. */
  GraphType graph_type =
      query_context_->graph_metadata->numPartitions() > 1
          ? GraphType::GraphView
          : (query_context_->query_config.use_auxiliary_index ? GraphType::BipartiteGraphView : GraphType::Normal);
  planner_ = std::make_unique<NaivePlanner>(&query_context_->query_graph, &cardinality, graph_type);

  // now order is directly configured
  auto use_order = getOrder(query_context_->query_config.matching_order, query_context_->query_graph.getNumVertices());

  ExecutionPlan* plan = nullptr;
  bool inputs_are_keys = true;
  // TODO(tatiana): consider make the first vertex in the matching order to be a key if multithread?
  switch (query_context_->query_config.compression_strategy) {
  case CompressionStrategy::Static: {
    plan = planner_->generatePlan(use_order);
    inputs_are_keys = plan->isInCover(plan->getRootQueryVertexID());
    break;
  }
  case CompressionStrategy::Dynamic: {
    planner_->generateOrder(use_order);
    std::vector<std::vector<double>> cardinality_per_level(cardinality.size(), cardinality);
    planner_->generateCoverNode(cardinality_per_level);
    plan = planner_->generatePlanWithDynamicCover();
    inputs_are_keys =
        plan->isInCover(plan->getRootQueryVertexID()) && (plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0);
    break;
  }
  case CompressionStrategy::None: {
    plan = planner_->generatePlanWithoutCompression(use_order);
    inputs_are_keys = true;
  }
  }

  LOG(INFO) << "Generated plan";
  backtracking_plan_ = std::make_unique<BacktrackingPlan>(
      plan, inputs_are_keys, plan->getRootQueryVertexID());  // now assume all vertices have candidates
  return backtracking_plan_.get();
}

}  // namespace circinus

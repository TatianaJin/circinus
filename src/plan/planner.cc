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
#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
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
}

CandidatePruningPlan* Planner::updateCandidatePruningPlan(const std::vector<VertexID>* cardinality) {
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
    switch (strategy) {
    case CandidatePruningStrategy::DAF: {
      break;
    }
    case CandidatePruningStrategy::CFL: {
      break;
    }
    case CandidatePruningStrategy::GQL: {
      break;
    }
    case CandidatePruningStrategy::TSO: {
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

ExecutionPlan* Planner::generateExecutionPlan(const std::vector<VertexID>*) {
  // TODO(tatiana)
  return nullptr;
}

}  // namespace circinus

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

#include <string>

#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "utils/query_utils.h"

namespace circinus {

CandidatePruningPlan* Planner::generateCandidatePruningPlan() {
  auto& q = query_context_->query_graph;
  auto& strategy = query_context_->query_config.candidate_pruning_strategy;
  if (strategy == "adaptive") {
    LOG(FATAL) << "Not implemented yet";
    // TODO(tatiana): adapt to data graph distribution
    return &candidate_pruning_plan_;
  }
  if (strategy == "ldf") {
    candidate_pruning_plan_.newLDFScan(q);
    return &candidate_pruning_plan_;
  }
  if (strategy == "dpiso") {
    candidate_pruning_plan_.newLDFScan(q);
    // FIXME(tatiana): add phase 2 filters for dpiso
    return &candidate_pruning_plan_;
  }
  if (strategy == "cfl") {
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
    // FIXME(tatiana): add phase 2 filters for cfl
    return &candidate_pruning_plan_;
  }
  if (strategy == "tso") {
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
    // FIXME(tatiana): add phase 2 filters for tso
    return &candidate_pruning_plan_;
  }
  if (strategy == "gql") {
    candidate_pruning_plan_.newLDFScan(q);
    candidate_pruning_plan_.newNLFFilter(q);
    // FIXME(tatiana): add phase 2 filters for tso
    return &candidate_pruning_plan_;
  }
  CHECK(false) << "Unknown strategy " << strategy;
  return &candidate_pruning_plan_;
}

ExecutionPlan* Planner::generateExecutionPlan() {
  // FIXME(tatiana)
}

}  // namespace circinus

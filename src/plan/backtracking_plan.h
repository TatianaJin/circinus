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

#include <memory>
#include <utility>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph_metadata.h"
#include "graph/types.h"
#include "ops/input_operator.h"
#include "ops/logical/compressed_input.h"
#include "plan/execution_plan.h"
#include "plan/operator_tree.h"

namespace circinus {

class BacktrackingPlan {
  std::vector<ExecutionPlan*> plans_;
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> partitioned_plans_;
  std::vector<std::unique_ptr<LogicalCompressedInputOperator>> input_operators_;  // size = plans_.size();

 public:
  OperatorTree& getOperatorTree(uint32_t plan_idx = 0) { return plans_[plan_idx]->getOperators(); }

  std::unique_ptr<InputOperator> getInputOperator(uint32_t plan_idx, const GraphMetadata& metadata,
                                                  ExecutionConfig& exec_conf) {
    return nullptr;  // TODO(byli)
  }

  inline const auto& getPlans() const { return plans_; }

  inline uint32_t getNumPartitionedPlans() const { return partitioned_plans_.size(); }
  inline const auto& getPartitionedPlan(uint32_t idx) const { return partitioned_plans_[idx]; }

  inline uint32_t addPlan(ExecutionPlan* plan) {
    plans_.push_back(plan);
    return plans_.size() - 1;
  }

  /**
   * @param partitioned_plans A vector of {plan_idx, candiate scope for each query vertex}.
   */
  inline void addPartitionedPlans(std::vector<std::pair<uint32_t, std::vector<CandidateScope>>>&& partitioned_plans) {
    partitioned_plans_ = std::move(partitioned_plans);
  }

  inline void addInputOperator(std::unique_ptr<LogicalCompressedInputOperator>&& op) {
    input_operators_.push_back(std::move(op));
  }

  // FIXME(tatiana): deprecate the following functions?
  inline bool inputsAreKeys(uint32_t plan_idx = 0) const { return plans_[plan_idx]->inputAreKeys(); }

  inline uint32_t getInputCandidateIndex(uint32_t plan_idx = 0) const {
    return plans_[plan_idx]->getRootQueryVertexID();  // now assume all vertices have candidates
  }

  inline Outputs& getOutputs(uint32_t plan_idx = 0) const { return plans_[plan_idx]->getOutputs(); }
};

}  // namespace circinus

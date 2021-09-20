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
#include <string>
#include <utility>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph_metadata.h"
#include "graph/types.h"
#include "ops/input_operator.h"
#include "ops/logical/compressed_input.h"
#include "plan/execution_plan.h"

namespace circinus {

class BacktrackingPlan {
  // logical
  std::vector<ExecutionPlan*> plans_;
  std::vector<std::unique_ptr<LogicalCompressedInputOperator>> input_operators_;  // size = plans_.size();
  // partitioned and parallel
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> partitioned_plans_;
  std::vector<bool> segment_mask_;

 public:
  std::vector<Operator*>& getOperators(uint32_t plan_idx = 0) { return plans_[plan_idx]->getOperators(); }

  inline Operator* getOutputOperator(uint32_t plan_idx = 0) { return plans_[plan_idx]->getOperators().back(); }

  std::unique_ptr<InputOperator> getInputOperator(uint32_t plan_idx = 0) {
    return input_operators_[plan_idx]->toPhysicalOperators();
  }

  inline const auto& getPlans() const { return plans_; }
  inline const auto& getPlan(uint32_t idx) const { return plans_[idx]; }
  inline auto getNumInputOperators() const { return input_operators_.size(); }

  inline uint32_t getNumPartitionedPlans() const { return partitioned_plans_.size(); }
  inline const auto& getPartitionedPlan(uint32_t idx) const { return partitioned_plans_[idx]; }

  inline bool toSegment(uint32_t idx) const {
    DCHECK_LT(idx, segment_mask_.size());
    return segment_mask_[idx];
  }

  // TODO(engineering): merge plans with the same compression and order for better log/profile readabiliity?
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

  inline void setPartitionedPlanSegmentMask(std::vector<bool>&& segment_mask) {
    segment_mask_ = std::move(segment_mask);
  }

  inline void addInputOperator(std::unique_ptr<LogicalCompressedInputOperator>&& op) {
    input_operators_.push_back(std::move(op));
  }

  inline void replaceInputOperator(uint32_t idx, std::unique_ptr<LogicalCompressedInputOperator>&& op) {
    input_operators_[idx] = std::move(op);
  }

  std::string toString() const {
    std::stringstream ss;
    for (uint32_t i = 0; i < plans_.size(); ++i) {
      ss << "[ Plan " << i << " ]\n";
      ss << "0," << input_operators_[i]->toPhysicalOperators()->toString() << std::endl;
      plans_[i]->printPhysicalPlan(ss);
    }
    if (!partitioned_plans_.empty()) {
      ss << "[ Partitions ]\n";
    }
    for (auto& plan : partitioned_plans_) {
      ss << "Plan " << plan.first;
      for (auto& scope : plan.second) {
        ss << ' ';
        scope.print(ss);
      }
      ss << '\n';
    }
    return ss.str();
  }
};

}  // namespace circinus

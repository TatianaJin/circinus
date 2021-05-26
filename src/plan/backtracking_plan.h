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

#include <vector>

#include "plan/execution_plan.h"
#include "plan/operator_tree.h"

namespace circinus {

class BacktrackingPlan {
  ExecutionPlan* plan_;
  bool inputs_are_keys_;
  uint32_t input_candidate_index_;

 public:
  BacktrackingPlan(ExecutionPlan* plan, bool inputs_are_keys, uint32_t input_candidate_index)
      : plan_(plan), inputs_are_keys_(inputs_are_keys), input_candidate_index_(input_candidate_index) {}
  OperatorTree& getOperatorTree() { return plan_->getOperators(); }

  bool inputsAreKeys() const { return inputs_are_keys_; }
  uint32_t getInputCandidateIndex() const { return input_candidate_index_; }
  Outputs& getOutputs() const { return plan_->getOutputs(); }
};

}  // namespace circinus

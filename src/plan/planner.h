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

#include <string>

#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "utils/query_utils.h"

namespace circinus {

class Planner {
  QueryContext* const query_context_ = nullptr;
  CandidatePruningPlan candidate_pruning_plan_;

 public:
  explicit Planner(QueryContext& query_context) : query_context_(&query_context) {}

  /**
   * A candidate pruning plan has three phases, which are all optional.
   * <li> Phase 1 uses a chain of local filters to generate candidate sets for a set of query vertices
   * independently.</li>
   * <li> Phase 2 expands the generated candidate sets to the candidates of other query vertices.</li>
   * <li> Phase 3 prunes generated candidate sets all together.</li>
   */
  CandidatePruningPlan* generateCandidatePruningPlan();

  ExecutionPlan* generateExecutionPlan();
};

}  // namespace circinus

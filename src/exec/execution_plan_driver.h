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

#include "exec/plan_driver.h"
#include "exec/result.h"
#include "exec/traverse_task.h"
#include "plan/backtracking_plan.h"
#include "utils/query_utils.h"

namespace circinus {

/** Supports partition-parallel execution.
 */
class ExecutionPlanDriver : public PlanDriver {
  BacktrackingPlan* plan_;
  ExecutionResult* result_;  // owned by ExecutorManager

 public:
  explicit ExecutionPlanDriver(BacktrackingPlan* plan) : plan_(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue, ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

/** Alternative to ExecutionPlanDriver, supports matching-parallel execution.
 */
class MatchingParallelExecutionPlanDriver : public PlanDriver {
  BacktrackingPlan* plan_;
  ExecutionResult* result_;  // owned by ExecutorManager

 public:
  explicit MatchingParallelExecutionPlanDriver(BacktrackingPlan* plan) : plan_(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue, ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

}  // namespace circinus

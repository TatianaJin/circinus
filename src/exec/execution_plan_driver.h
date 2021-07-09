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

#include "exec/plan_driver.h"
#include "exec/result.h"
#include "exec/task.h"
#include "exec/traverse_task.h"
#include "plan/backtracking_plan.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class ExecutionPlanDriverBase : public PlanDriver {
 protected:
  BacktrackingPlan* plan_;
  std::unique_ptr<CandidateResult> candidate_result_ = nullptr;
  ExecutionResult* result_;  // owned by ExecutorManager
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
  QueryType query_type_;

 public:
  explicit ExecutionPlanDriverBase(BacktrackingPlan* plan) : plan_(plan) {}
  virtual ~ExecutionPlanDriverBase() {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void finishPlan(ThreadsafeQueue<ServerEvent>* reply_queue) {
    result_->setElapsedExecutionTime(toSeconds(start_time_, std::chrono::high_resolution_clock::now()));
    result_->setCount();
    reply_queue->push(std::move(*finish_event_));
    reset();
  }

  inline void collectTaskInfo(TaskBase* task) const {
    result_->addEnumerateTime(task->getExecutionTime());
    result_->collect(task);
  }
};

/** Supports partition-parallel execution.
 */
class ExecutionPlanDriver : public ExecutionPlanDriverBase {
 public:
  explicit ExecutionPlanDriver(BacktrackingPlan* plan) : ExecutionPlanDriverBase(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue, ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

/** Alternative to ExecutionPlanDriver, supports matching-parallel execution.
 */
class MatchingParallelExecutionPlanDriver : public ExecutionPlanDriverBase {
 private:
  uint32_t batch_size_;
  std::unique_ptr<InputOperator> input_op_ = nullptr;
  std::vector<bool> task_depleted_;
  std::vector<CandidateSetView> candidates_;

 public:
  explicit MatchingParallelExecutionPlanDriver(BacktrackingPlan* plan) : ExecutionPlanDriverBase(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue, ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

}  // namespace circinus

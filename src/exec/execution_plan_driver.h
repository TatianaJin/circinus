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
#include <unordered_set>
#include <utility>
#include <vector>

#include "zmq.hpp"

#include "exec/executor_to_cost_learner_client.h"
#include "exec/plan_driver.h"
#include "exec/profile_task.h"
#include "exec/result.h"
#include "exec/task.h"
#include "exec/traverse_task.h"
#include "plan/backtracking_plan.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

#define DEFAULT_SUSPEND_INTERVAL 0.2

class ExecutionPlanDriverBase : public PlanDriver {
 protected:
  uint32_t batch_size_;
  BacktrackingPlan* plan_;
  std::unique_ptr<CandidateResult> candidate_result_ = nullptr;
  ExecutionResult* result_;  // owned by ExecutorManager
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
  QueryType query_type_;
  std::vector<std::vector<CandidateSetView>> candidates_;
  std::vector<std::vector<std::shared_ptr<unordered_set<VertexID>>>> candidate_hashmaps_;
  bool is_time_out_ = false;

  constexpr static float running_average_decay_ = 0.5;
  double running_average_enumerate_time_ = DEFAULT_SUSPEND_INTERVAL;
  double min_suspend_interval_ = DEFAULT_SUSPEND_INTERVAL;
  double suspend_interval_ = DEFAULT_SUSPEND_INTERVAL;  // determined by average task enumerate time and min interval
  double* suspend_interval_ptr_ = nullptr;
  uint32_t max_parallelism_ = 1;
  uint32_t n_finished_task_instances_ = 0;
  uint32_t n_pending_tasks_ = 0;

 public:
  explicit ExecutionPlanDriverBase(BacktrackingPlan* plan) : plan_(plan) {}
  virtual ~ExecutionPlanDriverBase() {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void finishPlan(ThreadsafeQueue<ServerEvent>* reply_queue);

  void taskTimeOut(std::unique_ptr<TaskBase>& task, ThreadsafeQueue<ServerEvent>* reply_queue) override;

 protected:
  inline void collectTaskInfo(std::unique_ptr<TaskBase>& task) {
    running_average_enumerate_time_ = running_average_enumerate_time_ * (1 - running_average_decay_) +
                                      running_average_decay_ * task->getExecutionTime();
    ++n_finished_task_instances_;
    result_->addEnumerateTime(task->getExecutionTime());
    result_->collect(task);
  }

  inline void addTaskToQueue(ThreadsafeTaskQueue* task_queue, std::unique_ptr<TaskBase>& task) {
    task_queue->putTask(std::move(task));
  }

  template <typename TaskType, typename... Args>
  inline void addTaskToQueue(ThreadsafeTaskQueue* task_queue, Args... args) {
    TaskBase* task = nullptr;
    if (isProfileMode(query_type_)) {
      task = new ProfileTask<TaskType>(std::forward<Args>(args)...);
    } else {
      task = new TaskType(std::forward<Args>(args)...);
    }
    if (std::is_same_v<TaskType, TraverseTask> || std::is_same_v<TaskType, TraverseChainTask>) {
      ((TraverseChainTask*)task)->setSuspendIntervalPtr(suspend_interval_ptr_);
      ((TraverseChainTask*)task)->setSplitSize(max_parallelism_);
    }
    task_queue->putTask(task);
  }

  template <typename TaskType>
  void handleSupendedTask(std::unique_ptr<TaskBase>&, ThreadsafeTaskQueue*);

  inline void updateSuspendInterval() {
    if (n_pending_tasks_ > max_parallelism_) {
      suspend_interval_ =
          std::max(min_suspend_interval_, running_average_enumerate_time_ * (n_pending_tasks_ / max_parallelism_));
    } else {
      suspend_interval_ = min_suspend_interval_;
    }
  }
};

/** Execution plan driver of online query with seed vertex
 */
class OnlineQueryExecutionPlanDriver : public ExecutionPlanDriverBase {
 protected:
 public:
  explicit OnlineQueryExecutionPlanDriver(BacktrackingPlan* plan) : ExecutionPlanDriverBase(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                  ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

/** Execution plan driver of online query with seed vertex
 */
class OnlineQueryExecutionPlanDriver : public ExecutionPlanDriverBase {
 public:
  OnlineQueryExecutionPlanDriver(BacktrackingPlan* plan) : ExecutionPlanDriverBase(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                  ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

/** Supports partition-parallel execution.
 */
class ExecutionPlanDriver : public ExecutionPlanDriverBase {
  QueryGraph* query_ = nullptr;
  uint32_t n_partitions_ = 0;
  uint32_t n_labels_ = 0;
  std::unique_ptr<ExecutorToCostLearnerClient> cost_learner_client_;

 public:
  ExecutionPlanDriver(BacktrackingPlan* plan, zmq::context_t* zmq_ctx) : ExecutionPlanDriverBase(plan) {
    if (!FLAGS_cost_learner.empty()) {
      CHECK_NOTNULL(zmq_ctx);
      cost_learner_client_ = std::make_unique<ExecutorToCostLearnerClient>(zmq_ctx, FLAGS_cost_learner);
    }
  }

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                  ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

/** Alternative to ExecutionPlanDriver, supports matching-parallel execution.
 */
class MatchingParallelExecutionPlanDriver : public ExecutionPlanDriverBase {
 private:
  std::unique_ptr<InputOperator> input_op_ = nullptr;
  std::vector<bool> task_depleted_;
  std::vector<CandidateSetView> candidates_;
  std::vector<std::unique_ptr<TraverseContext>> traverse_ctx_templates_;

 public:
  explicit MatchingParallelExecutionPlanDriver(BacktrackingPlan* plan) : ExecutionPlanDriverBase(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                  ThreadsafeQueue<ServerEvent>* reply_queue) override;

  void taskTimeOut(std::unique_ptr<TaskBase>& task, ThreadsafeQueue<ServerEvent>* reply_queue) override;
};

}  // namespace circinus

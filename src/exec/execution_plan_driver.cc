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

#include "exec/execution_plan_driver.h"

#include <memory>
#include <utility>

#include "exec/plan_driver.h"
#include "exec/result.h"
#include "exec/traverse_task.h"
#include "plan/execution_plan.h"
#include "utils/query_utils.h"

namespace circinus {

void ExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                               ThreadsafeTaskQueue& task_queue) {
  candidate_result_.reset(dynamic_cast<CandidateResult*>(ctx.second.release()));
  ctx.second = Result::newExecutionResult();
  result_ = (ExecutionResult*)ctx.second.get();
  result_->getOutputs().init(ctx.first.getNumExecutors()).limit(query_ctx->query_config.limit);

  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::ExecutionPhase);
  finish_event_->data = result_->data();
  finish_event_->query_id = qid;

  auto n_plans = plan_->getNumPartitionedPlans();
  CHECK_EQ(plan_->getPlans().size(), n_plans) << "now one partitioned plan per graph partition";
  task_counters_.resize(plan_->getPlans().size(), 1);  // one task per partition
  for (uint32_t i = 0; i < n_plans; ++i) {
    auto plan_idx = plan_->getPartitionedPlan(i).first;
    plan_->getOperatorTree(plan_idx).setOutput(&result_->getOutputs());  // all plans share the same output
    task_queue.putTask(new TraverseTask(qid, i, ctx.first.getBatchSize(), plan_->getOperatorTree(plan_idx),
                                        plan_->getInputOperator(plan_idx, *query_ctx->graph_metadata, ctx.first),
                                        plan_->getPartitionedPlan(i).second, query_ctx->data_graph,
                                        candidate_result_.get()));
  }
}

void ExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
    // TODO(tatiana): now we only consider count as output
    result_->setCount(result_->getOutputs().getCount());
    reply_queue->push(std::move(*finish_event_));
    reset();
  }
}

void MatchingParallelExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                               ThreadsafeTaskQueue& task_queue) {
  auto candidate_result = std::move(ctx.second);
  ctx.second = Result::newExecutionResult();
  result_ = (ExecutionResult*)ctx.second.get();
  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::ExecutionPhase);
  finish_event_->data = result_->data();
  finish_event_->query_id = qid;

  plan_->getOutputs().init(ctx.first.getNumExecutors()).limit(query_ctx->query_config.limit);
  if (ctx.first.getNumExecutors() == 1) {
    // one task for single-thread execution
    task_counters_.push_back(1);
    task_queue.putTask(new TraverseChainTask(
        qid, 0, ctx.first.getBatchSize(), plan_->getOperatorTree(), (const Graph*)query_ctx->data_graph,
        *dynamic_cast<CandidateResult*>(candidate_result.get())->getMergedCandidates(), plan_->getInputCandidateIndex(),
        plan_->inputsAreKeys()));
    return;
  }
  // // init task counter for root and put tasks to task queue
  // task_counters_.resize(plan_->getOperators().getOperatorSize());
  // auto root_ptr = plan_->getOperators().root();
  // task_counters_[0] = root_ptr->getParallelism();
  // for (uint32_t i = 0; i < root_ptr->getParallelism(); ++i) {
  //   auto traverse = dynamic_cast<TraverseOperator*>(root_ptr);
  //   task_queue.putTask(new TraverseTask(qid, 0, i, traverse, query_ctx->data_graph));
  // }
}

void MatchingParallelExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
    // TODO(tatiana): now we only consider count as output
    result_->setCount(plan_->getOutputs().getCount());
    reply_queue->push(std::move(*finish_event_));
    reset();
  }
  // TODO(tatiana)
}

}  // namespace circinus

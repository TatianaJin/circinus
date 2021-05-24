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

#include "exec/plan_driver.h"
#include "exec/result.h"
#include "exec/traverse_task.h"
#include "plan/execution_plan.h"
#include "utils/query_utils.h"

namespace circinus {

// TODO(tatiana): support partitioned graph
void ExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                               ThreadsafeTaskQueue& task_queue) {
  ctx.second = Result::newExecutionResult();
  result_ = (ExecutionResult*)ctx.second.get();
  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::ExecutionPhase);
  finish_event_->data = result_->data();
  finish_event_->query_id = qid;
  // init task counter for root and put tasks to task queue
  task_counters_.resize(plan_->getOperators().getOperatorSize());
  auto root_ptr = plan_->getOperators().root();
  task_counters_[0] = root_ptr->getParallelism();
  for (uint32_t i = 0; i < root_ptr->getParallelism(); ++i) {
    auto traverse = dynamic_cast<TraverseOperator*>(root_ptr);
    task_queue.putTask(new TraverseTask(qid, 0, i, traverse, query_ctx->data_graph));
  }
}

void ExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  // TODO(tatiana)
}

void MatchingParallelExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                               ThreadsafeTaskQueue& task_queue) {
  // TODO(tatiana)
}

void MatchingParallelExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  // TODO(tatiana)
}

}  // namespace circinus

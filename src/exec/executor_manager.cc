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

#include "exec/executor_manager.h"

#include <utility>

#include "exec/result.h"

namespace circinus {

ExecutorManager::ExecutorManager(ThreadsafeQueue<ServerEvent>* queue) : reply_queue_(queue) {
  executors_.start(&task_queue_, &finished_tasks_);
  finish_task_handler_ = std::thread([this]() {
    while (true) {
      auto task = finished_tasks_.waitAndPop();
      if (task == nullptr) return;
      PlanDriver* driver = nullptr;
      ExecutionContext* ctx = nullptr;
      {
        std::lock_guard<std::mutex> lock(execution_ctx_mu_);
        auto& state = execution_ctx_.at(task->getQueryId());
        ctx = &state.first;
        driver = state.second.get();
      }
      DCHECK(ctx->second != nullptr);
      ctx->second->collect(task.get());
      driver->taskFinish(task.get(), &task_queue_, reply_queue_);
    }
  });
}

void ExecutorManager::run(QueryId qid, QueryContext* query_ctx, std::unique_ptr<PlanDriver>&& plan_driver) {
  std::pair<ExecutionContext, std::unique_ptr<PlanDriver>>* ctx;
  {
    std::lock_guard<std::mutex> lock(execution_ctx_mu_);
    auto ctx_pos = execution_ctx_.find(qid);
    if (ctx_pos == execution_ctx_.end()) {
      ctx_pos =
          execution_ctx_
              .insert({qid, std::make_pair(std::make_pair(getExecutionConfig(), nullptr), std::move(plan_driver))})
              .first;
    } else if (plan_driver != nullptr) {
      ctx_pos->second.second = std::move(plan_driver);  // update driver
    }
    ctx = &ctx_pos->second;
  }
  ctx->second->init(qid, query_ctx, ctx->first, task_queue_);
}

void ExecutorManager::ExecutorPool::start(ThreadsafeTaskQueue* task_queue,
                                          ThreadsafeQueue<std::unique_ptr<TaskBase>>* finished_task) {
  pool_.reserve(n_threads_);
  for (uint32_t i = 0; i < n_threads_; ++i) {
    pool_.push_back(std::thread([i, task_queue, finished_task]() {
      while (true) {
        auto task = task_queue->getTask();
        task->run();
        finished_task->push(std::move(task));
      }
    }));
  }
}

// TODO(tatiana): setup ExecutionConfig
ExecutionConfig ExecutorManager::getExecutionConfig() const { return ExecutionConfig(); }

}  // namespace circinus

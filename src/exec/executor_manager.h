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
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "exec/plan_driver.h"
#include "exec/result.h"
#include "exec/task.h"
#include "exec/task_queue.h"
#include "exec/threadsafe_queue.h"
#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "utils/query_utils.h"

namespace circinus {

/** Accepts tasks from the circinus server and supports asynchronous execution of tasks.
 *
 * A thread pool is maintained to run tasks in parallel.
 */
class ExecutorManager {
 private:
  using ExecutionContext = std::pair<ExecutionConfig, std::unique_ptr<Result>>;

  class ExecutorPool {
   private:
    std::vector<std::thread> pool_;
    uint32_t n_threads_;

   public:
    explicit ExecutorPool(uint32_t n_threads = FLAGS_num_cores) : n_threads_(n_threads) {}

    ~ExecutorPool() {
      for (auto& t : pool_) {
        t.join();
      }
    }

    void start(ThreadsafeTaskQueue* task_queue, ThreadsafeQueue<std::unique_ptr<TaskBase>>* finished_task);
    inline void shutDown(ThreadsafeTaskQueue* task_queue) {
      for (uint32_t i = 0; i < n_threads_; ++i) {
        task_queue->putTask(nullptr);
      }
    }
    inline uint32_t getNumExecutors() const { return n_threads_; }
  } executors_;

  ThreadsafeQueue<ServerEvent>* reply_queue_;  // to server, not owned
  ThreadsafeTaskQueue task_queue_;
  ThreadsafeQueue<std::unique_ptr<TaskBase>> finished_tasks_;

  // stl is used as phmap invalidates the reference/pointer to element upon flat map mutation
  std::unordered_map<QueryId, std::pair<ExecutionContext, std::unique_ptr<PlanDriver>>> execution_ctx_;
  std::mutex execution_ctx_mu_;
  std::thread finish_task_handler_;

 public:
  explicit ExecutorManager(ThreadsafeQueue<ServerEvent>* queue);

  ~ExecutorManager() {
    executors_.shutDown(&task_queue_);
    finished_tasks_.push(nullptr);
    finish_task_handler_.join();
  }

  /* interface with circinus server, called by the circinus main thread */
  /**
   * @param plan_driver If plan_driver is nullptr, use the previous one.
   */
  void run(QueryId qid, QueryContext* query_ctx, std::unique_ptr<PlanDriver>&& plan_driver);
  inline void clearQuery(QueryId qid) {
    std::lock_guard<std::mutex> lock(execution_ctx_mu_);
    LOG(INFO) << "Query " << qid << " Finished.";
    execution_ctx_.erase(qid);
  }

 private:
  ExecutionConfig getExecutionConfig() const;
};

}  // namespace circinus

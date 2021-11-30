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

#include "exec/result.h"
#include "exec/task.h"
#include "exec/task_queue.h"
#include "exec/threadsafe_queue.h"
#include "plan/candidate_pruning_plan.h"
#include "utils/query_utils.h"

namespace circinus {

/** Plan Driver interface */
class PlanDriver {
 protected:
  /* for each query phase */
  std::vector<uint32_t> task_counters_;
  std::vector<std::unique_ptr<Operator>> operators_;
  TaskId n_finished_tasks_ = 0;
  std::unique_ptr<ServerEvent> finish_event_ = nullptr;

 public:
  using ExecutionContext = std::pair<ExecutionConfig, std::unique_ptr<Result>>;

  virtual ~PlanDriver() {
    for (auto& op : operators_) {
      LOG(INFO) << op->toString();
      google::FlushLogFiles(google::INFO);
      op.release();
    }
  }
  virtual void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) = 0;
  virtual void taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                          ThreadsafeQueue<ServerEvent>* reply_queue) = 0;
  virtual void taskTimeOut(std::unique_ptr<TaskBase>& task, ThreadsafeQueue<ServerEvent>* reply_queue) = 0;

 protected:
  void reset() {
    operators_.clear();
    task_counters_.clear();
    n_finished_tasks_ = 0;
  }
};

}  // namespace circinus

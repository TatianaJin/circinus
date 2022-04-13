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

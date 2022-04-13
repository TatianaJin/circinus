#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include "exec/filter_task.h"
#include "exec/plan_driver.h"
#include "exec/result.h"
#include "plan/candidate_pruning_plan.h"
#include "utils/query_utils.h"

namespace circinus {

class CandidatePruningPlanDriver : public PlanDriver {
  CandidatePruningPlan* plan_;
  CandidateResult* result_ = nullptr;  // owned by ExecutorManager
  ExecutionConfig* exec_config_ = nullptr;

 public:
  explicit CandidatePruningPlanDriver(CandidatePruningPlan* plan) : plan_(plan) {}

  void init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx, ThreadsafeTaskQueue& task_queue) override;

  void taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                  ThreadsafeQueue<ServerEvent>* reply_queue) override;

  void taskTimeOut(std::unique_ptr<TaskBase>& task, ThreadsafeQueue<ServerEvent>* reply_queue) override;

 private:
  void initPhase1TasksForPartitionedGraph(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                          ThreadsafeTaskQueue& task_queue);
  inline uint32_t getFilterParallelism(uint64_t input_size) const {
    auto batch_size = exec_config_->getBatchSize();
    uint32_t batches = (input_size + batch_size - 1) / batch_size;
    return std::min(batches, exec_config_->getNumExecutors());
  }
};

}  // namespace circinus

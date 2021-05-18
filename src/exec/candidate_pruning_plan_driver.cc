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

#include "exec/candidate_pruning_plan_driver.h"

#include <memory>
#include <vector>

#include "exec/filter_task.h"
#include "exec/plan_driver.h"
#include "exec/result.h"
#include "exec/scan_task.h"
#include "plan/candidate_pruning_plan.h"
#include "utils/query_utils.h"

namespace circinus {

// TODO(tatiana): support partitioned graph
void CandidatePruningPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                      ThreadsafeTaskQueue& task_queue) {
  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::CandidatePhase);
  finish_event_->data = &candidate_cardinality_;
  finish_event_->query_id = qid;
  if (plan_->getPhase() == 1) {
    auto scans = plan_->getScanOperators(*query_ctx->graph_metadata, ctx.first);
    task_counters_.resize(scans.size());
    ctx.second = Result::newCandidateResult(scans.size());
    result_ = (CandidateResult*)ctx.second.get();
    for (uint32_t task_id = 0; task_id < scans.size(); ++task_id) {
      auto& scan = scans[task_id];
      task_counters_[task_id] = scan->getParallelism();
      for (uint32_t i = 0; i < scan->getParallelism(); ++i) {
        task_queue.putTask(new ScanTask(qid, task_id, i, scan.get(), query_ctx->data_graph));
      }
      DLOG(INFO) << "Query " << qid << " Task " << task_id << " " << scan->getParallelism() << " shard(s)";
      operators_.push_back(std::move(scan));
    }
  } else if (plan_->getPhase() == 3) {
    auto filters = plan_->getFilterOperators(*query_ctx->graph_metadata, ctx.first);
    task_counters_.resize(filters.size());
    tasks_.resize(filters.size());
    for (uint32_t task_id = 0; task_id < filters.size(); ++task_id) {
      auto& filter = filters[task_id];
      task_counters_[task_id] = filter->getParallelism();
      tasks_[task_id].reserve(filter->getParallelism());
      for (uint32_t i = 0; i < filter->getParallelism(); ++i) {
        tasks_[task_id].emplace_back(std::make_unique<NeighborhoodFilterTask>(
            qid, task_id, i, filter.get(), query_ctx->data_graph, result_->getMergedCandidates()));
      }
    }
    for (uint32_t i = 0; i < tasks_[n_finished_tasks_].size(); ++i) {
      task_queue.putTask(tasks_[n_finished_tasks_][i].get());
    }
  }
}

void CandidatePruningPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                            ThreadsafeQueue<ServerEvent>* reply_queue) {
  if (--task_counters_[task->getTaskId()] == 0) {
    if (++n_finished_tasks_ == task_counters_.size()) {
      if (plan_->getPhase() == 1 || plan_->getPhase() == 2) {
        candidate_cardinality_ = result_->getCandidateCardinality();
        result_->merge();
        reply_queue->push(std::move(*finish_event_));
        reset();
      } else {
        result_->remove_invalid(tasks_[n_finished_tasks_ - 1][0]->getFilter()->getQueryVertex());
        reply_queue->push(std::move(*finish_event_));
        reset();
      }
      return;
    }

    if (plan_->getPhase() == 3) {
      result_->remove_invalid(tasks_[n_finished_tasks_ - 1][0]->getFilter()->getQueryVertex());
      for (uint32_t i = 0; i < tasks_[n_finished_tasks_].size(); ++i) {
        task_queue->putTask(tasks_[n_finished_tasks_][i].get());
      }
    }
  }
}

}  // namespace circinus

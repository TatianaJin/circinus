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
#include <utility>
#include <vector>

#include "exec/filter_task.h"
#include "exec/plan_driver.h"
#include "exec/result.h"
#include "exec/scan_task.h"
#include "plan/candidate_pruning_plan.h"
#include "utils/flags.h"
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
    uint32_t task_id = 0;
    auto& filter = filters[task_id];
    // Parallelism
    QueryVertexID query_vertex = filter->getQueryVertex();
    uint64_t input_size = (*result_->getMergedCandidates())[query_vertex].size();
    filter->setInputSize(input_size);
    filter->setParallelism((input_size + FLAGS_batch_size - 1) / FLAGS_batch_size);
    task_counters_[task_id] = filter->getParallelism();
    for (uint32_t i = 0; i < filter->getParallelism(); ++i) {
      task_queue.putTask(new NeighborhoodFilterTask(qid, task_id, i, filter.get(), query_ctx->data_graph,
                                                    result_->getMergedCandidates()));
    }
    for (auto& filter : filters) {
      operators_.push_back(std::move(filter));
    }
  }
}

void CandidatePruningPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                            ThreadsafeQueue<ServerEvent>* reply_queue) {
  // TODO(BYLI) package merge operation and remove_invalid operation into tasks, determine the parallelism
  if (--task_counters_[task->getTaskId()] == 0) {
    LOG(INFO) << "CandidatePhase " << plan_->getPhase() << " finished";
    if (++n_finished_tasks_ == task_counters_.size()) {
      if (plan_->getPhase() == 1 || plan_->getPhase() == 2) {
        result_->merge();
        candidate_cardinality_ = result_->getMergedCandidateCardinality();
        reply_queue->push(std::move(*finish_event_));
        reset();
      } else {
        result_->remove_invalid(dynamic_cast<NeighborhoodFilterTask*>(task)->getFilter()->getQueryVertex());
        candidate_cardinality_ = result_->getMergedCandidateCardinality();
        reply_queue->push(std::move(*finish_event_));
        reset();
      }
      return;
    }

    if (plan_->getPhase() == 3) {
      result_->remove_invalid(dynamic_cast<NeighborhoodFilterTask*>(task)->getFilter()->getQueryVertex());
      auto filter = dynamic_cast<NeighborhoodFilter*>(operators_[n_finished_tasks_].get());
      QueryVertexID query_vertex = filter->getQueryVertex();
      uint64_t input_size = (*result_->getMergedCandidates())[query_vertex].size();
      filter->setInputSize(input_size);
      filter->setParallelism((input_size + FLAGS_batch_size - 1) / FLAGS_batch_size);
      task_counters_[n_finished_tasks_] = filter->getParallelism();
      for (uint32_t i = 0; i < filter->getParallelism(); ++i) {
        task_queue->putTask(new NeighborhoodFilterTask(task->getQueryId(), n_finished_tasks_, i, filter,
                                                       task->getDataGraph(), result_->getMergedCandidates()));
      }
    }
  }
}

}  // namespace circinus

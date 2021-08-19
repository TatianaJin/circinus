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

void CandidatePruningPlanDriver::initPhase1TasksForPartitionedGraph(QueryId qid, QueryContext* query_ctx,
                                                                    ExecutionContext& ctx,
                                                                    ThreadsafeTaskQueue& task_queue) {
  auto& metadata = *query_ctx->graph_metadata;
  auto n_partitions = metadata.numPartitions();
  auto n_qvs = plan_->getScanQueryVertices().size();

  DCHECK_NE(n_qvs, 0) << "The case of no candidate is not handled yet";
  task_counters_.resize(n_qvs);
  ctx.second = Result::newPartitionedCandidateResult(n_qvs, n_partitions);
  result_ = (CandidateResult*)ctx.second.get();
  finish_event_->data = result_;
  // TODO(tatiana) now for simplicity enforce one task per partition for one query vertex
  ctx.first.setMaxParallelism(1);
  std::vector<ScanTask*> tasks;
  tasks.reserve(n_partitions * n_qvs);
  for (uint32_t i = 0; i < n_partitions; ++i) {
    auto scans = plan_->getScanOperators(metadata.getPartition(i), ctx.first);
    for (uint32_t task_id = 0; task_id < n_qvs; ++task_id) {
      auto& scan = scans[task_id];
      if (scan != nullptr) {
        task_counters_[task_id] += scan->getParallelism();
        tasks.push_back(new ScanTask(qid, task_id, query_ctx->stop_time, 0, scan.get(), query_ctx->data_graph, i));
        operators_.push_back(std::move(scan));
      } else {
        DLOG(INFO) << "partition " << i << '/' << n_partitions << " no candidates for " << task_id;
      }
      if (i == n_partitions - 1) {
        DLOG(INFO) << "Query " << qid << " Task " << task_id << " " << task_counters_[task_id] << " shard(s)";
      }
    }
  }
  for (auto counter : task_counters_) {
    if (counter == 0) {
      // TODO(tatiana): handle trivial case when there is no candidate
      LOG(FATAL) << " No candidate matching query vertex?";
    }
  }
  for (auto task : tasks) {
    task_queue.putTask(task);
  }
}

void CandidatePruningPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                      ThreadsafeTaskQueue& task_queue) {
  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::CandidatePhase);
  finish_event_->query_id = qid;
  if (plan_->getPhase() == 1) {
    LOG(INFO) << "Candidate pruning phase 1";
    auto n_partitions = query_ctx->graph_metadata->numPartitions();
    if (n_partitions > 1) {
      initPhase1TasksForPartitionedGraph(qid, query_ctx, ctx, task_queue);
      return;
    }
    auto scans = plan_->getScanOperators(*query_ctx->graph_metadata, ctx.first);
    DCHECK_NE(scans.size(), 0) << "The case of no candidate is not handled yet";
    task_counters_.resize(scans.size());
    ctx.second = Result::newCandidateResult(scans.size());
    result_ = (CandidateResult*)ctx.second.get();
    for (uint32_t task_id = 0; task_id < scans.size(); ++task_id) {
      auto& scan = scans[task_id];
      CHECK(scan != nullptr) << task_id;
      task_counters_[task_id] = scan->getParallelism();
      for (uint32_t i = 0; i < scan->getParallelism(); ++i) {
        task_queue.putTask(new ScanTask(qid, task_id, query_ctx->stop_time, i, scan.get(), query_ctx->data_graph));
      }
      DLOG(INFO) << "Query " << qid << " Task " << task_id << " " << scan->getParallelism() << " shard(s)";
      operators_.push_back(std::move(scan));
    }
  } else if (plan_->getPhase() == 3) {
    LOG(INFO) << "Candidate pruning phase 3";
    auto filters = plan_->getFilterOperators(*query_ctx->graph_metadata, ctx.first);
    CHECK_NE(filters.size(), 0) << "The case of no candidate is not handled yet";
    for (auto& filter : filters) {
      operators_.push_back(std::move(filter));
    }
    task_counters_.resize(filters.size());
    uint32_t task_id = 0;
    auto filter = dynamic_cast<NeighborhoodFilter*>(operators_[n_finished_tasks_].get());
    // Parallelism
    QueryVertexID query_vertex = filter->getQueryVertex();
    uint64_t input_size = (*result_->getMergedCandidates())[query_vertex].size();
    filter->setInputSize(input_size);
    filter->setParallelism((input_size + FLAGS_batch_size - 1) / FLAGS_batch_size);
    task_counters_[task_id] = filter->getParallelism();
    for (uint32_t i = 0; i < filter->getParallelism(); ++i) {
      task_queue.putTask(new NeighborhoodFilterTask(qid, task_id, query_ctx->stop_time, i, filter,
                                                    query_ctx->data_graph, result_->getMergedCandidates()));
    }
  }
  finish_event_->data = result_;
}

void CandidatePruningPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                            ThreadsafeQueue<ServerEvent>* reply_queue) {
  CHECK_LT(task->getTaskId(), task_counters_.size()) << "phase " << plan_->getPhase();
  if (plan_->getPhase() == 1) {
    // DLOG(INFO) << task->getTaskId() << ' ' << dynamic_cast<ScanTask*>(task)->getScanContext().candidates.size();
    result_->collect(task);
  }
  // TODO(BYLI) package merge operation and remove_invalid operation into tasks, determine the parallelism
  if (--task_counters_[task->getTaskId()] == 0) {
    if (plan_->getPhase() == 1) {
      if (!plan_->toPartitionResult()) {
        result_->merge(task);
      }
    }

    if (++n_finished_tasks_ == task_counters_.size()) {
      if (plan_->getPhase() == 3) {
        result_->removeInvalid(dynamic_cast<NeighborhoodFilterTask*>(task)->getFilter()->getQueryVertex());
      }
      reset();
      reply_queue->push(std::move(*finish_event_));  // must be the last step in phase to avoid data race
      return;
    }

    if (plan_->getPhase() == 3) {
      result_->removeInvalid(dynamic_cast<NeighborhoodFilterTask*>(task)->getFilter()->getQueryVertex());
      auto filter = dynamic_cast<NeighborhoodFilter*>(operators_[n_finished_tasks_].get());
      QueryVertexID query_vertex = filter->getQueryVertex();
      uint64_t input_size = (*result_->getMergedCandidates())[query_vertex].size();
      CHECK_NE(input_size, 0) << "query vertex " << query_vertex << " has empty candidate set, not handled yet";
      filter->setInputSize(input_size);
      filter->setParallelism((input_size + FLAGS_batch_size - 1) / FLAGS_batch_size);
      task_counters_[n_finished_tasks_] = filter->getParallelism();
      for (uint32_t i = 0; i < filter->getParallelism(); ++i) {
        task_queue->putTask(new NeighborhoodFilterTask(task->getQueryId(), n_finished_tasks_, task->getStopTime(), i,
                                                       filter, task->getDataGraph(), result_->getMergedCandidates()));
      }
    }
  }
}

void CandidatePruningPlanDriver::taskTimeOut(TaskBase* task, ThreadsafeQueue<ServerEvent>* reply_queue) {
  if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
    finish_event_->type = ServerEvent::TimeOut;
    reply_queue->push(std::move(*finish_event_));
  }
}

}  // namespace circinus

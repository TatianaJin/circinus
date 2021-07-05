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
#include <vector>

#include "exec/plan_driver.h"
#include "exec/profile_task.h"
#include "exec/result.h"
#include "exec/traverse_task.h"
#include "plan/execution_plan.h"
#include "utils/query_utils.h"

namespace circinus {

void ExecutionPlanDriverBase::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                   ThreadsafeTaskQueue& task_queue) {
  start_time_ = std::chrono::high_resolution_clock::now();
  candidate_result_.reset(dynamic_cast<CandidateResult*>(ctx.second.release()));
  ctx.second = Result::newExecutionResult(query_ctx->query_config.mode == QueryMode::Profile, plan_->getPlans().size());
  result_ = (ExecutionResult*)ctx.second.get();
  result_->getOutputs().init(ctx.first.getNumExecutors()).limit(query_ctx->query_config.limit);
  if (plan_->getPlans().size() == 1) {
    std::stringstream ss;
    for (auto qv : plan_->getPlans().front()->getMatchingOrder()) {
      ss << qv << ' ';
    }
    result_->setMatchingOrder(ss.str());
  } else if (plan_->getPlans().size() > 1) {
    result_->setMatchingOrder("mixed");
  }

  // TODO(tatiana): add profile mode ProfileWithMiniIntersection
  query_type_ = (query_ctx->query_config.mode == QueryMode::Profile)
                    ? QueryType::Profile
                    : (/*(query_ctx->query_config.mode == QueryMode::ProfileWithMiniIntersection) ?
                          QueryType::ProfileWithMiniIntersection :*/
                       QueryType::Execute);

  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::ExecutionPhase);
  finish_event_->data = &result_->getQueryResult();
  finish_event_->query_id = qid;
}

void ExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                               ThreadsafeTaskQueue& task_queue) {
  ExecutionPlanDriverBase::init(qid, query_ctx, ctx, task_queue);

  auto n_plans = plan_->getNumPartitionedPlans();
  CHECK_EQ(plan_->getPlans().size(), n_plans) << "now one partitioned plan per graph partition";
  task_counters_.resize(plan_->getPlans().size(), 1);  // one task per partition
  for (uint32_t i = 0; i < n_plans; ++i) {
    auto plan_idx = plan_->getPartitionedPlan(i).first;
    dynamic_cast<OutputOperator*>(plan_->getOutputOperator(plan_idx))
        ->setOutput(&result_->getOutputs());  // all plans share the same output
    auto input_operator = plan_->getInputOperator(plan_idx);

    auto& scopes = plan_->getPartitionedPlan(i).second;
    auto partitioned_result = dynamic_cast<PartitionedCandidateResult*>(candidate_result_.get());
    if (query_ctx->query_config.mode == QueryMode::Profile) {
      task_queue.putTask(new ProfileTask<TraverseTask>(qid, i, ctx.first.getBatchSize(), plan_->getOperators(plan_idx),
                                                       std::move(input_operator), scopes, query_ctx->data_graph,
                                                       partitioned_result->getCandidatesByScopes(scopes)));
    } else {
      task_queue.putTask(new TraverseTask(qid, i, ctx.first.getBatchSize(), plan_->getOperators(plan_idx),
                                          std::move(input_operator), scopes, query_ctx->data_graph,
                                          partitioned_result->getCandidatesByScopes(scopes)));
    }
  }
}

void ExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  collectTaskInfo(task);
  if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
    // TODO(tatiana): now we only consider count as output
    finishPlan(reply_queue);
  }
}

void MatchingParallelExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                               ThreadsafeTaskQueue& task_queue) {
  ExecutionPlanDriverBase::init(qid, query_ctx, ctx, task_queue);
  dynamic_cast<OutputOperator*>(plan_->getOutputOperator())->setOutput(&result_->getOutputs());

  result_->getOutputs().init(ctx.first.getNumExecutors()).limit(query_ctx->query_config.limit);
  if (false && ctx.first.getNumExecutors() == 1) {
    // one task for single-thread execution
    task_counters_.push_back(1);

    if (query_ctx->query_config.mode == QueryMode::Profile) {
      task_queue.putTask(new ProfileTask<TraverseChainTask>(qid, 0, ctx.first.getBatchSize(), plan_->getOperators(),
                                                            plan_->getInputOperator(), query_ctx->data_graph,
                                                            candidate_result_->getCandidates()));
    } else {
      task_queue.putTask(new TraverseChainTask(qid, 0, ctx.first.getBatchSize(), plan_->getOperators(),
                                               plan_->getInputOperator(), query_ctx->data_graph,
                                               candidate_result_->getCandidates()));
    }
    return;
  }

  {  // TODO(tatiana): wrap with a task?
    candidates_ = candidate_result_->getCandidates();

    for (auto op : plan_->getOperators()) {
      auto traverse = dynamic_cast<TraverseOperator*>(op);
      if (traverse == nullptr) break;
      traverse->setCandidateSets(&candidates_[traverse->getTargetQueryVertex()]);
    }
  }

  batch_size_ = ctx.first.getBatchSize();
  task_counters_.resize(plan_->getOperators().size() + 1, 0);  // input + traverse + output operator chain
  task_depleted_.resize(plan_->getOperators().size() + 1, false);
  input_op_ = plan_->getInputOperator();
  task_counters_[0] = input_op_->getParallelism();
  CHECK_EQ(task_counters_[0], 1) << " input task counter is not equal to 1 ";  // TODO(tatiana): parallel input tasks

  input_op_->setNext(plan_->getOperators().front());

  task_queue.putTask(new MatchingParallelInputTask(qid, 0, (const GraphBase*)query_ctx->data_graph, input_op_.get(),
                                                   candidate_result_->getCandidates()));
}

void MatchingParallelExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  collectTaskInfo(task);
  auto matching_parallel_task = dynamic_cast<MatchingParallelTask*>(task);
  if (matching_parallel_task == nullptr) {  // single-thread execution
    finishPlan(reply_queue);
    return;
  }

  auto op = matching_parallel_task->getNextOperator();
  if (op != nullptr) {
    auto inputs = std::make_shared<std::vector<CompressedSubgraphs>>(std::move(matching_parallel_task->getOutputs()));
    auto input_size = inputs->size();
    uint32_t level = matching_parallel_task->getNextLevel();
    task_counters_[level] += input_size / batch_size_ + (input_size % batch_size_ != 0);

    // Output Task
    if (task->getTaskId() == plan_->getOperators().size() - 1) {
      uint32_t input_index = 0;
      for (; input_index + batch_size_ < input_size; input_index += batch_size_) {
        task_queue->putTask(new MatchingParallelOutputTask(task->getQueryId(), level, dynamic_cast<OutputOperator*>(op),
                                                           inputs, input_index, input_index + batch_size_));
      }
      if (input_index < input_size) {
        task_queue->putTask(new MatchingParallelOutputTask(task->getQueryId(), level, dynamic_cast<OutputOperator*>(op),
                                                           inputs, input_index, input_size));
      }
    } else {
      uint32_t input_index = 0;
      for (; input_index + batch_size_ < input_size; input_index += batch_size_) {
        task_queue->putTask(new MatchingParallelTraverseTask(
            task->getQueryId(), level, dynamic_cast<TraverseOperator*>(op),
            dynamic_cast<TraverseOperator*>(op)->initTraverseContext(inputs.get(), task->getDataGraph(), input_index,
                                                                     input_index + batch_size_, query_type_),
            inputs));
      }
      if (input_index < input_size) {
        task_queue->putTask(new MatchingParallelTraverseTask(
            task->getQueryId(), level, dynamic_cast<TraverseOperator*>(op),
            dynamic_cast<TraverseOperator*>(op)->initTraverseContext(inputs.get(), task->getDataGraph(), input_index,
                                                                     input_size, query_type_),
            inputs));
      }
    }
  }

  if (--task_counters_[task->getTaskId()] == 0) {
    if (task->getTaskId() == 0) {
      task_depleted_[task->getTaskId()] = true;
    } else {
      task_depleted_[task->getTaskId()] = task_depleted_[task->getTaskId() - 1];
      // TODO(tatiana) support match limit
    }
    if (task_depleted_[task->getTaskId()] && ++n_finished_tasks_ == task_counters_.size()) {
      finishPlan(reply_queue);
    }
  }
}

}  // namespace circinus

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

#include <chrono>
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
  ctx.second = Result::newExecutionResult(query_ctx->query_config.isProfileMode(), plan_->getPlans());
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

  query_type_ = (query_ctx->query_config.mode == QueryMode::Profile)
                    ? QueryType::Profile
                    : ((query_ctx->query_config.mode == QueryMode::ProfileSI) ? QueryType::ProfileWithMiniIntersection
                                                                              : QueryType::Execute);
  LOG(INFO) << "Query Type " << ((uint16_t)query_type_);

  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::ExecutionPhase);
  finish_event_->data = &result_->getQueryResult();
  finish_event_->query_id = qid;
}

void ExecutionPlanDriverBase::finishPlan(ThreadsafeQueue<ServerEvent>* reply_queue) {
  if (is_time_out_) {
    finish_event_->type = ServerEvent::TimeOut;
    reply_queue->push(std::move(*finish_event_));
    reset();
    return;
  }

  result_->setElapsedExecutionTime(toSeconds(start_time_, std::chrono::high_resolution_clock::now()));
  result_->setCount();
  auto profiled = dynamic_cast<ProfiledExecutionResult*>(result_);
  if (profiled != nullptr) {
    auto size = plan_->getPlans().size();
    std::stringstream ss;
    for (uint32_t i = 0; i < size; ++i) {
      ss << "[ Plan " << plan_->getPlan(i)->getPartitionId() << " ]\n";
      auto input_op = plan_->getInputOperator(i);
      profiled->setProfiledPlan(i, plan_->getOperators(i), input_op.get());
      uint32_t op_idx = 0;
      for (auto& op_profile : profiled->getProfiledPlanStrings()[i]) {
        ss << op_idx << ',' << op_profile << std::endl;
        ++op_idx;
      }
    }
    if (plan_->getNumPartitionedPlans() != 0) {
      ss << "[ Partitions ]\n";
    }
    for (uint32_t i = 0; i < plan_->getNumPartitionedPlans(); ++i) {
      auto& plan = plan_->getPartitionedPlan(i);
      ss << "Plan " << plan.first;
      for (auto& scope : plan.second) {
        ss << ' ';
        scope.print(ss);
      }
      ss << '\n';
    }
    result_->getQueryResult().profile_info = ss.str();
  }
  reply_queue->push(std::move(*finish_event_));
  reset();
}

void ExecutionPlanDriverBase::taskTimeOut(TaskBase* task, ThreadsafeQueue<ServerEvent>* reply_queue) {
  is_time_out_ = true;
  if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
    finishPlan(reply_queue);
  }
}

void ExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                               ThreadsafeTaskQueue& task_queue) {
  ExecutionPlanDriverBase::init(qid, query_ctx, ctx, task_queue);

  auto n_plans = plan_->getNumPartitionedPlans();
  // FIXME(tatiana): rank partitioned plans by their weights
  task_counters_.resize(n_plans, 1);  // one task per partition
  candidates_.resize(n_plans);
  for (uint32_t i = 0; i < n_plans; ++i) {
    // if (i != 5) continue;
    auto plan_idx = plan_->getPartitionedPlan(i).first;
    // all plans share the same output
    dynamic_cast<OutputOperator*>(plan_->getOutputOperator(plan_idx))->setOutput(&result_->getOutputs());
    // plan_->getPlan(plan_idx)->printPhysicalPlan();
    auto input_operator = plan_->getInputOperator(plan_idx);

    auto& scopes = plan_->getPartitionedPlan(i).second;
    auto partitioned_result = dynamic_cast<PartitionedCandidateResult*>(candidate_result_.get());
    candidates_[i] = partitioned_result->getCandidatesByScopes(scopes);
    addTaskToQueue<TraverseTask>(&task_queue, qid, i, query_ctx->stop_time, ctx.first.getBatchSize(),
                                 plan_->getOperators(plan_idx), std::move(input_operator), scopes,
                                 query_ctx->data_graph, &candidates_[i], query_type_);
  }
}

void ExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  collectTaskInfo(task);
  if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
    finishPlan(reply_queue);
  }
}

void MatchingParallelExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                               ThreadsafeTaskQueue& task_queue) {
  ExecutionPlanDriverBase::init(qid, query_ctx, ctx, task_queue);
  dynamic_cast<OutputOperator*>(plan_->getOutputOperator())->setOutput(&result_->getOutputs());

  result_->getOutputs().init(ctx.first.getNumExecutors()).limit(query_ctx->query_config.limit);

  candidates_ = candidate_result_->getCandidates();
  if (ctx.first.getNumExecutors() == 1) {
    // one task for single-thread execution
    task_counters_.push_back(1);

    addTaskToQueue<TraverseChainTask>(&task_queue, qid, 0, query_ctx->stop_time, ctx.first.getBatchSize(),
                                      plan_->getOperators(), plan_->getInputOperator(), query_ctx->data_graph,
                                      &candidates_, query_type_);
    return;
  }

  {  // TODO(tatiana): wrap with a task?
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
  addTaskToQueue<MatchingParallelInputTask>(&task_queue, qid, 0, query_ctx->stop_time,
                                            (const GraphBase*)query_ctx->data_graph, input_op_.get(),
                                            candidate_result_->getCandidates());
}

void MatchingParallelExecutionPlanDriver::taskFinish(TaskBase* task, ThreadsafeTaskQueue* task_queue,
                                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  collectTaskInfo(task);
  auto matching_parallel_task = dynamic_cast<MatchingParallelTask*>(task);
  if (matching_parallel_task == nullptr) {  // single-thread execution
    finishPlan(reply_queue);
    return;
  }

  // LOG(INFO) << "Task " << task->getTaskId() << " Finish.";
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
        addTaskToQueue<MatchingParallelOutputTask>(task_queue, task->getQueryId(), level, task->getStopTime(),
                                                   dynamic_cast<OutputOperator*>(op), inputs, input_index,
                                                   input_index + batch_size_);
      }
      if (input_index < input_size) {
        addTaskToQueue<MatchingParallelOutputTask>(task_queue, task->getQueryId(), level, task->getStopTime(),
                                                   dynamic_cast<OutputOperator*>(op), inputs, input_index, input_size);
      }
    } else {
      uint32_t input_index = 0;
      for (; input_index + batch_size_ < input_size; input_index += batch_size_) {
        addTaskToQueue<MatchingParallelTraverseTask>(
            task_queue, task->getQueryId(), level, task->getStopTime(), dynamic_cast<TraverseOperator*>(op),
            dynamic_cast<TraverseOperator*>(op)->initTraverseContext(inputs.get(), task->getDataGraph(), input_index,
                                                                     input_index + batch_size_, query_type_),
            inputs);
      }
      if (input_index < input_size) {
        addTaskToQueue<MatchingParallelTraverseTask>(
            task_queue, task->getQueryId(), level, task->getStopTime(), dynamic_cast<TraverseOperator*>(op),
            dynamic_cast<TraverseOperator*>(op)->initTraverseContext(inputs.get(), task->getDataGraph(), input_index,
                                                                     input_size, query_type_),
            inputs);
      }
    }
  }

  if (--task_counters_[task->getTaskId()] == 0) {
    // TODO(limit): abort task when match limit is reached
    if (task->getTaskId() == 0 || task_depleted_[task->getTaskId() - 1]) {
      for (uint32_t i = task->getTaskId(); i < task_counters_.size() && task_counters_[i] == 0; ++i) {
        task_depleted_[i] = true;
        ++n_finished_tasks_;
      }
    }
    if (n_finished_tasks_ == task_counters_.size()) {
      finishPlan(reply_queue);
    }
  }
}

void MatchingParallelExecutionPlanDriver::taskTimeOut(TaskBase* task, ThreadsafeQueue<ServerEvent>* reply_queue) {
  is_time_out_ = true;
  if (--task_counters_[task->getTaskId()] == 0) {
    if (task->getTaskId() == 0 || task_depleted_[task->getTaskId() - 1]) {
      for (uint32_t i = task->getTaskId(); task_counters_[i] == 0 && i < task_counters_.size(); ++i) {
        task_depleted_[i] = true;
        ++n_finished_tasks_;
      }
    }
    if (n_finished_tasks_ == task_counters_.size()) {
      finishPlan(reply_queue);
    }
  }
}

}  // namespace circinus

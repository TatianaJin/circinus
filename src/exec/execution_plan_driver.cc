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

#include <algorithm>
#include <chrono>
#include <memory>
#include <unordered_set>
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
  batch_size_ = ctx.first.getBatchSize();
  max_parallelism_ = ctx.first.getNumExecutors();
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

  switch (query_ctx->query_config.mode) {
  case QueryMode::Profile:
    query_type_ = QueryType::Profile;
    break;
  case QueryMode::ProfileSI:
    query_type_ = QueryType::ProfileWithMiniIntersection;
    break;
  case QueryMode::ProfileCandidateSI:
    query_type_ = QueryType::ProfileCandidateSIEffect;
    break;
  default:
    query_type_ = QueryType::Execute;
  }
  if (verboseExecutionLog()) {
    LOG(INFO) << "Query Type " << ((uint16_t)query_type_);
  }

  finish_event_ = std::make_unique<ServerEvent>(ServerEvent::ExecutionPhase);
  finish_event_->data = &result_->getQueryResult();
  finish_event_->query_id = qid;
}

void ExecutionPlanDriverBase::finishPlan(ThreadsafeQueue<ServerEvent>* reply_queue) {
  if (is_time_out_) {
    finish_event_->type = ServerEvent::TimeOut;
    reset();
    reply_queue->push(std::move(*finish_event_));  // must be the last step in phase to avoid data race
    return;
  }

  LOG(INFO) << "num finished task instances " << n_finished_task_instances_ << " running_average_enumerate_time_ "
            << running_average_enumerate_time_;
  result_->setElapsedExecutionTime(toSeconds(start_time_, std::chrono::high_resolution_clock::now()));
  result_->setCount();
  auto profiled = dynamic_cast<ProfiledExecutionResult*>(result_);
  if (profiled != nullptr) {
    auto size = plan_->getPlans().size();
    std::stringstream ss;
    for (uint32_t i = 0; i < size; ++i) {
      ss << "[ Plan " << i << ": " << (plan_->getOperators(i).size() + 1) << " ]\n";
      auto input_op = plan_->getInputOperator(i);
      profiled->setProfiledPlan(i, plan_->getOperators(i), input_op.get());
      uint32_t op_idx = 0;
      for (auto& op_profile : profiled->getProfiledPlanStrings()[i]) {
        ss << op_idx << ',' << op_profile << std::endl;
        ++op_idx;
      }
      ss << "step_costs";
      for (auto step_cost : plan_->getPlan(i)->getStepCosts()) {
        ss << ' ' << step_cost;
      }
      ss << std::endl;
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
  reset();
  reply_queue->push(std::move(*finish_event_));  // must be the last step in phase to avoid data race
}

void ExecutionPlanDriverBase::taskTimeOut(std::unique_ptr<TaskBase>& task, ThreadsafeQueue<ServerEvent>* reply_queue) {
  is_time_out_ = true;
  if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
    finishPlan(reply_queue);
  }
}

template <typename TaskType>
void ExecutionPlanDriverBase::handleSupendedTask(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue) {
  auto traverse_task = dynamic_cast<TaskType*>(task.get());
  auto splits = traverse_task->getSplits();
  CHECK(!splits.empty());
  n_pending_tasks_ += splits.size();
  task_counters_[task->getTaskId()] += splits.size();
  uint32_t qid = traverse_task->getQueryId();
  uint32_t task_id = traverse_task->getTaskId();
  std::chrono::time_point<std::chrono::steady_clock> stop_time = traverse_task->getStopTime();
  auto& ops = plan_->getOperators(task_id);
  if (std::is_same_v<TaskType, TraverseTask>) {
    for (auto& split : splits) {
      auto input_size = split.second.size();
      addTaskToQueue<TraverseTask>(task_queue, qid, task_id, stop_time, batch_size_, &ops,
                                   ((TraverseTask*)traverse_task)->getScopes(), traverse_task->getDataGraph(),
                                   traverse_task->getCandidates(), query_type_, split.first, ops.size() - 1,
                                   std::move(split.second), 0, input_size, traverse_task->getCandidateHashmaps());
    }
  } else {
    for (auto& split : splits) {
      auto input_size = split.second.size();
      addTaskToQueue<TraverseChainTask>(task_queue, qid, task_id, stop_time, batch_size_, &ops,
                                        traverse_task->getDataGraph(), traverse_task->getCandidates(), query_type_,
                                        split.first, ops.size() - 1, std::move(split.second), 0, input_size,
                                        traverse_task->getCandidateHashmaps());
    }
  }
  if (traverse_task->getTaskStatus() == TaskStatus::Suspended) {
    // put into the queue for continuing execution
    addTaskToQueue(task_queue, task);
  } else {
    --n_pending_tasks_;
    --task_counters_[task->getTaskId()];
  }
  updateSuspendInterval();
}

void OnlineQueryExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                                          ThreadsafeTaskQueue& task_queue) {
  ExecutionPlanDriverBase::init(qid, query_ctx, ctx, task_queue);
  if (query_ctx->query_config.seed.first != DUMMY_QUERY_VERTEX) {
    VertexID seed_data_vertex = query_ctx->query_config.seed.second;
    task_counters_.resize(1, 1);
    n_pending_tasks_ += 1;
    dynamic_cast<OutputOperator*>(plan_->getOutputOperator(0))->setOutput(&result_->getOutputs());
    const std::vector<Operator*>& ops = plan_->getOperators(0);
    addTaskToQueue<TraverseChainTask>(&task_queue, qid, 0, query_ctx->stop_time, batch_size_, &ops,
                                      query_ctx->data_graph, query_type_, seed_data_vertex);
  } else {
    auto n_plans = plan_->getNumPartitionedPlans();
    task_counters_.resize(plan_->getPlans().size(), 0);
    candidates_.resize(n_plans);
    for (uint32_t i = 0; i < n_plans; ++i) {
      ++task_counters_[plan_->getPartitionedPlan(i).first];
    }
    n_pending_tasks_ += n_plans;

    for (uint32_t i = 0; i < n_plans; ++i) {
      auto plan_idx = plan_->getPartitionedPlan(i).first;
      // all plans share the same output
      dynamic_cast<OutputOperator*>(plan_->getOutputOperator(plan_idx))->setOutput(&result_->getOutputs());
      std::unique_ptr<InputOperator> input_operator = std::move(plan_->getInputOperator(plan_idx));
      auto& scopes = plan_->getPartitionedPlan(i).second;
      CHECK_EQ(scopes.size(), 1) << "Expecting scope only for input (start) query vertex";
      const std::vector<Operator*>& ops = plan_->getOperators(plan_idx);
      auto partitioned_result = dynamic_cast<PartitionedCandidateResult*>(candidate_result_.get());
      DCHECK(partitioned_result != nullptr);
      candidates_[i] = partitioned_result->getCandidatesByScopes(scopes);

      if (max_parallelism_ > 1 && plan_->toSegment(i)) {
        suspend_interval_ptr_ = &suspend_interval_;
      } else {
        suspend_interval_ptr_ = nullptr;
      }

      addTaskToQueue<TraverseChainTask>(&task_queue, qid, plan_idx, query_ctx->stop_time, batch_size_, &ops,
                                        std::move(input_operator), query_ctx->data_graph, &candidates_[i], query_type_,
                                        ops.size() - 1, nullptr);
    }
  }
}

void OnlineQueryExecutionPlanDriver::taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                                                ThreadsafeQueue<ServerEvent>* reply_queue) {
  auto traverse_task = dynamic_cast<TraverseChainTask*>(task.get());
  if (traverse_task->getTaskStatus() != TaskStatus::Normal) {
    handleSupendedTask<TraverseChainTask>(task, task_queue);
  } else {
    --n_pending_tasks_;
    collectTaskInfo(task);
    updateSuspendInterval();
    if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
      DCHECK_EQ(n_pending_tasks_, 0);
      finishPlan(reply_queue);
    }
  }
}

void ExecutionPlanDriver::init(QueryId qid, QueryContext* query_ctx, ExecutionContext& ctx,
                               ThreadsafeTaskQueue& task_queue) {
  ExecutionPlanDriverBase::init(qid, query_ctx, ctx, task_queue);

  // for cost learner
  query_ = &query_ctx->query_graph;
  n_partitions_ = query_ctx->query_config.use_partitioned_graph ? query_ctx->graph_metadata->numPartitions() : 1;
  n_labels_ = query_ctx->data_graph->getNumLabels();

  auto n_plans = plan_->getNumPartitionedPlans();
  // TODO(tatiana): handle the case with no plan (no match after pruning)
  CHECK_NE(n_plans, 0) << "The case of no plan is not handled yet";
  task_counters_.resize(plan_->getPlans().size(), 0);
  candidates_.resize(n_plans);
  for (uint32_t i = 0; i < n_plans; ++i) {
    ++task_counters_[plan_->getPartitionedPlan(i).first];
  }
  n_pending_tasks_ += n_plans;

  // >>>>> hard code now to share hashmaps of candidate sets
  candidate_hashmaps_.resize(n_plans);
  // FIXME(tatiana): problematic if using LDF/NLF cps as candidates are not merged
  auto candidate_views = candidate_result_->getCandidates();
  CHECK_EQ(candidate_views.size(), query_ctx->query_graph.getNumVertices())
      << "Now assuming generating candidates for all qvs";
  std::vector<std::shared_ptr<unordered_set<VertexID>>> full_candidate_hashmaps;
  full_candidate_hashmaps.reserve(candidate_views.size());
  for (auto& view : candidate_views) {
    full_candidate_hashmaps.push_back(std::make_shared<unordered_set<VertexID>>(view.begin(), view.end()));
  }
  // <<<<< hard code now to share hashmaps of candidate sets

  for (uint32_t i = 0; i < n_plans; ++i) {
    auto plan_idx = plan_->getPartitionedPlan(i).first;
    // all plans share the same output
    dynamic_cast<OutputOperator*>(plan_->getOutputOperator(plan_idx))->setOutput(&result_->getOutputs());
    // plan_->getPlan(plan_idx)->printPhysicalPlan();
    std::unique_ptr<InputOperator> input_operator = std::move(plan_->getInputOperator(plan_idx));

    auto& scopes = plan_->getPartitionedPlan(i).second;
    const std::vector<Operator*>& ops = plan_->getOperators(plan_idx);

    auto partitioned_result = dynamic_cast<PartitionedCandidateResult*>(candidate_result_.get());
    DCHECK(partitioned_result != nullptr);
    candidates_[i] = partitioned_result->getCandidatesByScopes(scopes);

    // >>>>> hard code now to share hashmaps of candidate sets
    candidate_hashmaps_[i].resize(scopes.size(), nullptr);
    for (uint32_t qv = 0; qv < scopes.size(); ++qv) {
      if (scopes[qv].getType() != CandidateScopeType::All) {
        auto& set = candidates_[i][qv];
        candidate_hashmaps_[i][qv] = std::make_shared<unordered_set<VertexID>>(set.begin(), set.end());
      } else {
        candidate_hashmaps_[i][qv] = full_candidate_hashmaps[qv];
      }
    }
    // <<<<< hard code now to share hashmaps of candidate sets

    if (max_parallelism_ > 1 && plan_->toSegment(i)) {
      suspend_interval_ptr_ = &suspend_interval_;
    } else {
      suspend_interval_ptr_ = nullptr;
    }

    uint32_t end_level = ops.size() - 1;
    addTaskToQueue<TraverseTask>(&task_queue, qid, plan_idx, query_ctx->stop_time, ctx.first.getBatchSize(), &ops,
                                 std::move(input_operator), scopes, query_ctx->data_graph, &candidates_[i], query_type_,
                                 end_level, &candidate_hashmaps_[i]);
  }
}

void ExecutionPlanDriver::taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  DCHECK_LT(task->getTaskId(), task_counters_.size());
  auto traverse_task = dynamic_cast<TraverseTask*>(task.get());
  if (traverse_task->getTaskStatus() != TaskStatus::Normal) {
    handleSupendedTask<TraverseTask>(task, task_queue);
  } else {
    --n_pending_tasks_;
    collectTaskInfo(task);
    updateSuspendInterval();
    if (--task_counters_[task->getTaskId()] == 0 && ++n_finished_tasks_ == task_counters_.size()) {
      if (cost_learner_client_ != nullptr) {
        auto profiled = dynamic_cast<ProfiledExecutionResult*>(result_);
        if (profiled != nullptr) {
          cost_learner_client_->sendUpdates(profiled, plan_, query_, n_partitions_, n_labels_, candidates_);
        }
      }
      DCHECK_EQ(n_pending_tasks_, 0);
      finishPlan(reply_queue);
    }
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
                                      &plan_->getOperators(), plan_->getInputOperator(), query_ctx->data_graph,
                                      &candidates_, query_type_);
    return;
  }

  {  // TODO(tatiana): wrap with a task?
    traverse_ctx_templates_.reserve(plan_->getOperators().size() - 1);
    for (auto op : plan_->getOperators()) {
      auto traverse = dynamic_cast<TraverseOperator*>(op);
      if (traverse == nullptr) break;
      traverse_ctx_templates_.emplace_back(
          traverse->initTraverseContext(&candidates_[traverse->getTargetQueryVertex()], nullptr,
                                        (const GraphBase*)query_ctx->data_graph, query_type_));
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

void MatchingParallelExecutionPlanDriver::taskFinish(std::unique_ptr<TaskBase>& task, ThreadsafeTaskQueue* task_queue,
                                                     ThreadsafeQueue<ServerEvent>* reply_queue) {
  collectTaskInfo(task);
  auto matching_parallel_task = dynamic_cast<MatchingParallelTask*>(task.get());
  if (matching_parallel_task == nullptr) {  // single-thread execution
    finishPlan(reply_queue);
    return;
  }

  LOG(INFO) << "Task " << task->getTaskId() << " Finish.";
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
        auto ctx_clone = traverse_ctx_templates_[level - 1]->clone();
        ctx_clone->setInput(*inputs, input_index, input_index + batch_size_);
        addTaskToQueue<MatchingParallelTraverseTask>(task_queue, task->getQueryId(), level, task->getStopTime(),
                                                     dynamic_cast<TraverseOperator*>(op), std::move(ctx_clone), inputs);
      }
      if (input_index < input_size) {
        auto ctx_clone = traverse_ctx_templates_[level - 1]->clone();
        ctx_clone->setInput(*inputs, input_index, input_size);
        addTaskToQueue<MatchingParallelTraverseTask>(task_queue, task->getQueryId(), level, task->getStopTime(),
                                                     dynamic_cast<TraverseOperator*>(op), std::move(ctx_clone), inputs);
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

void MatchingParallelExecutionPlanDriver::taskTimeOut(std::unique_ptr<TaskBase>& task,
                                                      ThreadsafeQueue<ServerEvent>* reply_queue) {
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

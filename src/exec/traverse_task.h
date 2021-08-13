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

#include <stdint.h>
#include <memory>
#include <utility>
#include <vector>

#include "glog/logging.h"
#ifdef WITH_GPERF
#include "gperftools/profiler.h"
#endif

#include "exec/task.h"
#include "graph/graph.h"
#include "graph/graph_partition.h"
#include "graph/graph_view.h"
#include "ops/input_operator.h"
#include "ops/output_operator.h"
#include "ops/traverse_operator.h"
#include "ops/types.h"
#include "utils/flags.h"

namespace circinus {

/**
 * For single-thread query execution
 */
class TraverseChainTask : public TaskBase {
 protected:
  const GraphBase* graph_ = nullptr;
  const uint32_t batch_size_;
  std::vector<Operator*> operators_;
  std::unique_ptr<InputOperator> input_op_ = nullptr;
  const std::vector<CandidateSetView>* candidates_;

  std::vector<ProfileInfo> profile_info_;
  QueryType query_type_;
  std::vector<std::unique_ptr<TraverseContext>> traverse_context_;
  std::vector<std::vector<CompressedSubgraphs>> outputs_;

 public:
  TraverseChainTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                    uint32_t batch_size, const std::vector<Operator*>& ops, std::unique_ptr<InputOperator>&& input_op,
                    const GraphBase* graph, const std::vector<CandidateSetView>* candidates, QueryType query_type)
      : TaskBase(qid, tid, stop_time),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        input_op_(std::move(input_op)),
        candidates_(candidates),
        query_type_(query_type) {
    CHECK(!candidates_->empty());
  }

  const GraphBase* getDataGraph() const override { return graph_; }
  const auto& getProfileInfo() const { return profile_info_; }

  void run(uint32_t executor_idx) override {
    auto start = std::chrono::steady_clock::now();
    setupOutputs();
    setupTraverseContexts();
    if (verboseExecutionLog()) {
      std::stringstream ss;
      for (auto& x : *candidates_) {
        ss << x.size() << " ";
      }
      LOG(INFO) << "-------- Task " << task_id_ << " candidate sizes " << ss.str();
    }
    auto old_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    // TODO(tatiana): support match limit
    auto inputs = input_op_->getInputs(graph_, *candidates_);
    // auto profile_output = "Task_" + std::to_string(task_id_);
    // ProfilerStart(profile_output.data());
    execute<QueryType::Execute>(inputs, inputs.size(), 0, executor_idx);
    // ProfilerStop();
    auto new_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    auto end = std::chrono::steady_clock::now();
    if (shortExecutionLog()) {
      LOG(INFO) << "Task " << task_id_ << " input " << inputs.size() << '/' << getNumSubgraphs(inputs, 0, inputs.size())
                << " count " << (new_count - old_count) << " (" << old_count << " to " << new_count
                << "), time usage: " << toSeconds(start, end) << "s.";
    }
  }

  void profile(uint32_t executor_idx) override {
    setupOutputs();
    setupTraverseContexts();
    LOG(INFO) << "Task " << task_id_;
    profile_info_.resize(operators_.size() + 1);  // traverse chain + input operator
    std::vector<CompressedSubgraphs> inputs;

    input_op_->inputAndProfile(graph_, *candidates_, &inputs, &profile_info_.front());
    LOG(INFO) << "Task " << task_id_;
    // auto profile_output = "Profile_Task_" + std::to_string(task_id_);
    // ProfilerStart(profile_output.data());
    // TODO(tatiana): support match limit
    auto start = std::chrono::steady_clock::now();
    auto old_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    LOG(INFO) << "Task " << task_id_;
    execute<QueryType::Profile>(inputs, inputs.size(), 0, executor_idx);
    for (uint32_t i = 0; i + 1 < operators_.size(); ++i) {
      profile_info_[i + 1] += *traverse_context_[i];
    }
    auto new_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    auto end = std::chrono::steady_clock::now();
    LOG(INFO) << "Task " << task_id_ << " input " << inputs.size() << " count " << (new_count - old_count)
              << ", time usage: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0
              << "s.";
    // ProfilerStop();
  }

 protected:
  void setupTraverseContexts() {
    auto op = operators_[0];
    auto traverse = dynamic_cast<TraverseOperator*>(op);
    uint32_t level = 0;
    traverse_context_.reserve(operators_.size() - 1);
    while (traverse != nullptr) {
      DCHECK_LT(traverse->getTargetQueryVertex(), candidates_->size());
      // now assume all query vertices have candidate sets
      traverse_context_.emplace_back(createTraverseContext(&(*candidates_)[traverse->getTargetQueryVertex()],
                                                           outputs_[level], level, traverse, query_type_));
      // LOG(INFO) << "set candidates for " << traverse->getTargetQueryVertex() << " " << traverse->toString();
      traverse = dynamic_cast<TraverseOperator*>(traverse->getNext());
      ++level;
    }
  }

  void setupOutputs() {
    auto output_size = input_op_->getOutputSize();
    outputs_.resize(operators_.size() - 1);
    for (uint32_t i = 0; i < outputs_.size(); ++i) {
      output_size = ((const TraverseOperator*)operators_[i])->getOutputSize(output_size);
      outputs_[i].resize(batch_size_, CompressedSubgraphs(output_size.first, output_size.second));
    }
  }

  virtual std::unique_ptr<TraverseContext> createTraverseContext(const CandidateSetView* candidates,
                                                                 std::vector<CompressedSubgraphs>& outputs,
                                                                 uint32_t level, const TraverseOperator* op,
                                                                 QueryType profile) {
    auto ctx = op->initTraverseContext(candidates, &outputs, graph_, profile);
    return ctx;
  }

  template <QueryType profile>
  bool execute(const std::vector<CompressedSubgraphs>& inputs, uint32_t input_size, uint32_t level,
               uint32_t executor_idx) {
    if (isTimeOut()) {
      return true;
    }
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      if
        constexpr(isProfileMode(profile)) {
          return output_op->validateAndOutputAndProfile(inputs, 0, input_size, executor_idx, &profile_info_[level + 1]);
        }
      uint32_t start = 0;
      return output_op->validateAndOutput(inputs, start, input_size, executor_idx);
    }

    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    auto ctx = traverse_context_[level].get();
    DCHECK_NOTNULL(ctx->getOutputs());
    ctx->setInput(inputs, 0, input_size);
    bool finished = false;
    while (true) {
      if (isTimeOut()) {
        return true;
      }
      ctx->resetOutputs();
      uint32_t size = 0;
      if
        constexpr(isProfileMode(profile)) size = traverse_op->expandAndProfile(batch_size_, ctx);
      else
        size = traverse_op->expand(batch_size_, ctx);

      if (size == 0) {
        break;
      }
      if (execute<profile>(outputs_[level], size, level + 1, executor_idx)) {
        finished = true;
        break;
      }
    }
    if
      constexpr(isProfileMode(profile)) { ctx->total_input_size += ctx->getTotalInputSize(); }

    return finished;
  }
};

class TraverseTask : public TraverseChainTask {
 private:
  const ReorderedPartitionedGraph* partitioned_graph_;
  std::vector<GraphView<GraphPartitionBase>> graph_views_;

 public:
  TraverseTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
               uint32_t batch_size, const std::vector<Operator*>& operators, std::unique_ptr<InputOperator>&& input_op,
               const std::vector<CandidateScope>& scopes, const GraphBase* graph,
               const std::vector<CandidateSetView>* candidates, QueryType query_type)
      : TraverseChainTask(query_id, task_id, stop_time, batch_size, operators, std::move(input_op), graph, candidates,
                          query_type),
        partitioned_graph_(dynamic_cast<const ReorderedPartitionedGraph*>(graph)) {
    CHECK(partitioned_graph_ != nullptr);
    graph_views_ = setupGraphView(partitioned_graph_, operators, scopes);
  }

 protected:
  std::unique_ptr<TraverseContext> createTraverseContext(const CandidateSetView* candidates,
                                                         std::vector<CompressedSubgraphs>& outputs, uint32_t level,
                                                         const TraverseOperator* op, QueryType profile) override {
    auto ctx = op->initTraverseContext(candidates, &outputs, &graph_views_[level], profile);
    return ctx;
  }

  static std::vector<GraphView<GraphPartitionBase>> setupGraphView(const ReorderedPartitionedGraph* g,
                                                                   const std::vector<Operator*>& operators,
                                                                   const std::vector<CandidateScope>& scopes) {
    std::vector<GraphView<GraphPartitionBase>> data_graphs_for_operators;
    CHECK_GT(operators.size(), 0);
    size_t len = operators.size() - 1;
    data_graphs_for_operators.reserve(len);
    for (size_t i = 0; i < len; ++i) {
      auto op = operators[i];
      auto traverse_op = dynamic_cast<TraverseOperator*>(op);
      CHECK(traverse_op != nullptr);
      data_graphs_for_operators.emplace_back(traverse_op->computeGraphPartitions(g, scopes));
    }
    return data_graphs_for_operators;
  }
};

/**
 * For multi-thread query execution
 */
class MatchingParallelTask : public TaskBase {
 protected:
  std::vector<CompressedSubgraphs> outputs_;

 public:
  MatchingParallelTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time)
      : TaskBase(qid, tid, stop_time) {}

  inline uint32_t getNextLevel() const { return getTaskId() + 1; }
  inline std::vector<CompressedSubgraphs>& getOutputs() { return outputs_; }

  virtual Operator* getNextOperator() = 0;
  virtual void collectProfileInfo(ProfileInfo& agg) const {}
};

class MatchingParallelInputTask : public MatchingParallelTask {
  const GraphBase* graph_ = nullptr;
  const InputOperator* input_op_ = nullptr;
  std::vector<CandidateSetView> candidates_;
  ProfileInfo info_;

 public:
  MatchingParallelInputTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                            const GraphBase* graph, const InputOperator* input_op,
                            std::vector<CandidateSetView>&& candidates)
      : MatchingParallelTask(qid, tid, stop_time),
        graph_(graph),
        input_op_(input_op),
        candidates_(std::move(candidates)) {}

  const GraphBase* getDataGraph() const override { return graph_; }

  void run(uint32_t executor_idx) override { outputs_ = input_op_->getInputs(graph_, candidates_); }

  void profile(uint32_t executor_idx) override { input_op_->inputAndProfile(graph_, candidates_, &outputs_, &info_); }

  void collectProfileInfo(ProfileInfo& agg) const override { agg += info_; }

  Operator* getNextOperator() override { return input_op_->getNext(); }
};

/**
 * For multi-thread query execution
 */
class MatchingParallelTraverseTask : public MatchingParallelTask {
  const TraverseOperator* traverse_op_;
  std::unique_ptr<TraverseContext> traverse_ctx_;
  std::vector<CandidateSetView> candidates_;
  std::shared_ptr<std::vector<CompressedSubgraphs>> inputs_;
  uint32_t batch_size_ = FLAGS_batch_size;  // FIXME(tatiana): use execution config

 public:
  MatchingParallelTraverseTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                               TraverseOperator* traverse_op, std::unique_ptr<TraverseContext>&& traverse_ctx,
                               const std::shared_ptr<std::vector<CompressedSubgraphs>>& inputs)
      : MatchingParallelTask(qid, tid, stop_time),
        traverse_op_(std::move(traverse_op)),
        traverse_ctx_(std::move(traverse_ctx)),
        inputs_(inputs) {
    DCHECK(traverse_op_ != nullptr);
    traverse_ctx_->setOutputBuffer(outputs_);
  }

  const GraphBase* getDataGraph() const override { return traverse_ctx_->getDataGraph<GraphBase>(); }

  void collectProfileInfo(ProfileInfo& agg) const override {
    traverse_ctx_->total_input_size = traverse_ctx_->getTotalInputSize();
    agg += *traverse_ctx_;
  }

  void run(uint32_t executor_idx) override {
    uint64_t output_size = 0;
    auto[nkey, size] = traverse_op_->getOutputSize(
        {traverse_ctx_->getCurrentInput().getNumKeys(), traverse_ctx_->getCurrentInput().getNumVertices()});
    CompressedSubgraphs placeholder(nkey, size);
    while (true) {
      outputs_.resize((output_size + batch_size_) / batch_size_ * batch_size_, placeholder);
      auto add = traverse_op_->expand(batch_size_, traverse_ctx_.get());
      if (add == 0) break;
      output_size += add;
    }
    outputs_.erase(outputs_.begin() + traverse_ctx_->getOutputSize(), outputs_.end());
  }

  void profile(uint32_t executor_idx) override {
    uint64_t output_size = 0;
    auto[nkey, size] = traverse_op_->getOutputSize(
        {traverse_ctx_->getCurrentInput().getNumKeys(), traverse_ctx_->getCurrentInput().getNumVertices()});
    CompressedSubgraphs placeholder(nkey, size);
    while (true) {
      outputs_.resize((output_size + batch_size_) / batch_size_ * batch_size_, placeholder);
      auto add = traverse_op_->expandAndProfile(batch_size_, traverse_ctx_.get());
      if (add == 0) break;
      output_size += add;
    }
    outputs_.erase(outputs_.begin() + traverse_ctx_->getOutputSize(), outputs_.end());
    traverse_ctx_->total_input_size = traverse_ctx_->getTotalInputSize();
  }

  Operator* getNextOperator() override { return traverse_op_->getNext(); }
};

class MatchingParallelOutputTask : public MatchingParallelTask {
  const OutputOperator* output_op_;
  ProfileInfo info_;
  std::shared_ptr<std::vector<CompressedSubgraphs>> inputs_;
  uint32_t input_index_ = 0;
  uint32_t input_end_index_ = 0;

 public:
  MatchingParallelOutputTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                             OutputOperator* op, const std::shared_ptr<std::vector<CompressedSubgraphs>>& inputs,
                             uint32_t input_index, uint32_t input_end_index)
      : MatchingParallelTask(qid, tid, stop_time),
        output_op_(op),
        inputs_(inputs),
        input_index_(input_index),
        input_end_index_(input_end_index) {
    DCHECK(output_op_ != nullptr);
  }

  const GraphBase* getDataGraph() const override { return nullptr; }

  void collectProfileInfo(ProfileInfo& agg) const override { agg += info_; }

  void run(uint32_t executor_idx) override {
    output_op_->validateAndOutput(*inputs_, input_index_, input_end_index_, executor_idx);
  }

  void profile(uint32_t executor_idx) override {
    output_op_->validateAndOutputAndProfile(*inputs_, input_index_, input_end_index_, executor_idx, &info_);
  }

  Operator* getNextOperator() override { return nullptr; }
};

}  // namespace circinus

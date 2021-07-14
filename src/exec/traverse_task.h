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

 public:
  TraverseChainTask(QueryId qid, TaskId tid, uint32_t batch_size, const std::vector<Operator*>& ops,
                    std::unique_ptr<InputOperator>&& input_op, const GraphBase* graph,
                    const std::vector<CandidateSetView>* candidates, QueryType query_type)
      : TaskBase(qid, tid),
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
    // auto profile_output = "Task " + std::to_string(task_id_);
    // ProfilerStart(profile_output.data());
    setupCandidateSets();

    auto old_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    // TODO(tatiana): support match limit
    auto inputs = input_op_->getInputs(graph_, *candidates_);
    execute<QueryType::Execute>(inputs, 0, executor_idx);
    auto new_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    LOG(INFO) << "Task " << task_id_ << " input " << inputs.size() << " count " << (new_count - old_count);
    // ProfilerStop();
  }

  void profile(uint32_t executor_idx) override {
    setupCandidateSets();
    profile_info_.resize(operators_.size() + 1);  // traverse chain + input operator
    std::vector<CompressedSubgraphs> inputs;
    input_op_->inputAndProfile(graph_, *candidates_, &inputs, &profile_info_.front());
    // TODO(tatiana): support match limit
    execute<QueryType::Profile>(inputs, 0, executor_idx);
  }

 protected:
  void setupCandidateSets() const {
    auto op = operators_[0];
    auto traverse = dynamic_cast<TraverseOperator*>(op);
    while (traverse != nullptr) {
      DCHECK_LT(traverse->getTargetQueryVertex(), candidates_->size());
      // now assume all query vertices have candidate sets
      traverse->setCandidateSets(&(*candidates_)[traverse->getTargetQueryVertex()]);
      // LOG(INFO) << "set candidates for " << traverse->getTargetQueryVertex() << " " << traverse->toString();
      op = op->getNext();
      traverse = dynamic_cast<TraverseOperator*>(op);
    }
  }

  virtual std::unique_ptr<TraverseContext> createTraverseContext(const std::vector<CompressedSubgraphs>& inputs,
                                                                 std::vector<CompressedSubgraphs>& outputs,
                                                                 uint32_t level, const TraverseOperator* op,
                                                                 QueryType profile) {
    auto ctx = op->initTraverseContext(&inputs, graph_, 0, inputs.size(), profile);
    ctx->outputs = &outputs;
    return ctx;
  }

  template <QueryType profile>
  bool execute(const std::vector<CompressedSubgraphs>& inputs, uint32_t level, uint32_t executor_idx) {
    std::vector<CompressedSubgraphs> outputs;
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      if
        constexpr(isProfileMode(profile)) {
          return output_op->validateAndOutputAndProfile(inputs, 0, inputs.size(), executor_idx,
                                                        &profile_info_[level + 1]);
        }
      return output_op->validateAndOutput(inputs, executor_idx);
    }
    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    auto ctx = createTraverseContext(inputs, outputs, level, traverse_op, query_type_);
    bool finished = false;
    while (true) {
      outputs.clear();
      uint32_t size = 0;
      if
        constexpr(isProfileMode(profile)) size = traverse_op->expandAndProfile(batch_size_, ctx.get());
      else
        size = traverse_op->expand(batch_size_, ctx.get());

      if (size == 0) {
        break;
      }
      if (execute<profile>(outputs, level + 1, executor_idx)) {
        finished = true;
        break;
      }
    }
    if
      constexpr(isProfileMode(profile)) {
        ctx->total_input_size = ctx->getTotalInputSize();
        profile_info_[level + 1] += *ctx;
      }

    return finished;
  }
};

class TraverseTask : public TraverseChainTask {
 private:
  const ReorderedPartitionedGraph* partitioned_graph_;
  std::vector<GraphView<GraphPartitionBase>> graph_views_;

 public:
  TraverseTask(QueryId query_id, TaskId task_id, uint32_t batch_size, const std::vector<Operator*>& operators,
               std::unique_ptr<InputOperator>&& input_op, const std::vector<CandidateScope>& scopes,
               const GraphBase* graph, const std::vector<CandidateSetView>* candidates, QueryType query_type)
      : TraverseChainTask(query_id, task_id, batch_size, operators, std::move(input_op), graph, candidates, query_type),
        partitioned_graph_(dynamic_cast<const ReorderedPartitionedGraph*>(graph)) {
    CHECK(partitioned_graph_ != nullptr);
    graph_views_ = setupGraphView(partitioned_graph_, operators, scopes);
  }

 protected:
  std::unique_ptr<TraverseContext> createTraverseContext(const std::vector<CompressedSubgraphs>& inputs,
                                                         std::vector<CompressedSubgraphs>& outputs, uint32_t level,
                                                         const TraverseOperator* op, QueryType profile) override {
    auto ctx = op->initTraverseContext(&inputs, &graph_views_[level], 0, inputs.size(), profile);
    ctx->outputs = &outputs;
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
  MatchingParallelTask(QueryId qid, TaskId tid) : TaskBase(qid, tid) {}

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
  MatchingParallelInputTask(QueryId qid, TaskId tid, const GraphBase* graph, const InputOperator* input_op,
                            std::vector<CandidateSetView>&& candidates)
      : MatchingParallelTask(qid, tid), graph_(graph), input_op_(input_op), candidates_(std::move(candidates)) {}

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

 public:
  MatchingParallelTraverseTask(QueryId qid, TaskId tid, TraverseOperator* traverse_op,
                               std::unique_ptr<TraverseContext>&& traverse_ctx,
                               const std::shared_ptr<std::vector<CompressedSubgraphs>>& inputs)
      : MatchingParallelTask(qid, tid),
        traverse_op_(std::move(traverse_op)),
        traverse_ctx_(std::move(traverse_ctx)),
        inputs_(inputs) {
    DCHECK(traverse_op_ != nullptr);
    traverse_ctx_->outputs = &outputs_;
  }

  const GraphBase* getDataGraph() const override {
    return reinterpret_cast<const GraphBase*>(traverse_ctx_->current_data_graph);
  }

  void collectProfileInfo(ProfileInfo& agg) const override {
    traverse_ctx_->total_input_size = traverse_ctx_->getTotalInputSize();
    agg += *traverse_ctx_;
  }

  void run(uint32_t executor_idx) override { traverse_op_->expand(UINT32_MAX, traverse_ctx_.get()); }
  void profile(uint32_t executor_idx) override {
    traverse_op_->expandAndProfile(UINT32_MAX, traverse_ctx_.get());
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
  MatchingParallelOutputTask(QueryId qid, TaskId tid, OutputOperator* op,
                             const std::shared_ptr<std::vector<CompressedSubgraphs>>& inputs, uint32_t input_index,
                             uint32_t input_end_index)
      : MatchingParallelTask(qid, tid),
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

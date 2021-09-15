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
#include <algorithm>
#include <memory>
#include <unordered_set>
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

enum class TaskStatus : uint64_t { Normal, Suspended, Split };
/**
 * For single-thread query execution
 */
class TraverseChainTask : public TaskBase {
 protected:
  using CandidateHashmaps = std::vector<std::shared_ptr<unordered_set<VertexID>>>;

  const GraphBase* graph_ = nullptr;
  const uint32_t batch_size_;
  const std::vector<Operator*>* operators_;
  const std::unique_ptr<InputOperator> input_op_ = nullptr;
  const std::vector<CandidateSetView>* candidates_ = nullptr;
  const CandidateHashmaps* candidate_hashmaps_ = nullptr;
  QueryType query_type_;
  uint32_t start_level_ = 0;
  uint32_t end_level_;
  TaskStatus task_status_ = TaskStatus::Normal;
  std::vector<CompressedSubgraphs> inputs_;
  uint32_t start_index_ = 0;
  uint32_t input_size_ = 0;

  // states for traverse
  std::vector<ProfileInfo> profile_info_;
  std::vector<std::unique_ptr<TraverseContext>> traverse_context_;
  std::vector<std::vector<CompressedSubgraphs>> outputs_;

  // dynamic segment
  uint32_t suspended_level_ = 0;
  double* suspend_interval_ = nullptr;
  uint32_t split_level_ = 0;
  uint32_t split_size_ = 1;
  uint32_t split_batch_size_ = 32;
  std::vector<std::pair<uint32_t, std::vector<CompressedSubgraphs>>> splits_;

  VertexID seed_data_vertex_ = INVALID_VERTEX_ID;

 public:
  // for online task with seed data vertex
  TraverseChainTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                    uint32_t batch_size, const std::vector<Operator*>* ops, const GraphBase* graph,
                    QueryType query_type, VertexID seed_data_vertex)
      : TaskBase(qid, tid, stop_time),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        query_type_(query_type),
        end_level_(operators_->size() - 1),
        inputs_(std::vector<CompressedSubgraphs>{CompressedSubgraphs(seed_data_vertex)}),
        input_size_(1) {}

  // for online task with seed data vertex
  TraverseChainTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                    uint32_t batch_size, const std::vector<Operator*>& ops, const GraphBase* graph,
                    QueryType query_type, VertexID seed_data_vertex)
      : TaskBase(qid, tid, stop_time),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        query_type_(query_type),
        start_level_(0),
        end_level_(operators_.size() - 1),
        task_status_(TaskStatus::Normal),
        seed_data_vertex_(seed_data_vertex) {}

  // for head task
  TraverseChainTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                    uint32_t batch_size, const std::vector<Operator*>* ops, std::unique_ptr<InputOperator>&& input_op,
                    const GraphBase* graph, const std::vector<CandidateSetView>* candidates, QueryType query_type,
                    uint32_t end_level, const CandidateHashmaps* hashmaps)
      : TaskBase(qid, tid, stop_time),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        input_op_(std::move(input_op)),
        candidates_(candidates),
        candidate_hashmaps_(hashmaps),
        query_type_(query_type),
        end_level_(end_level) {
    CHECK(!candidates_->empty());
  }

  // for head task
  TraverseChainTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                    uint32_t batch_size, const std::vector<Operator*>* ops, std::unique_ptr<InputOperator>&& input_op,
                    const GraphBase* graph, const std::vector<CandidateSetView>* candidates, QueryType query_type)
      : TraverseChainTask(qid, tid, stop_time, batch_size, ops, std::move(input_op), graph, candidates, query_type,
                          ops->size() - 1, nullptr) {
    CHECK(!candidates_->empty());
  }

  // for task starting from partial embeddings
  TraverseChainTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                    uint32_t batch_size, const std::vector<Operator*>* ops, const GraphBase* graph,
                    const std::vector<CandidateSetView>* candidates, QueryType query_type, uint32_t start_level,
                    uint32_t end_level, std::vector<CompressedSubgraphs>&& inputs, uint32_t start_index,
                    uint32_t input_size, const CandidateHashmaps* hashmaps)
      : TaskBase(qid, tid, stop_time),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        candidates_(candidates),
        candidate_hashmaps_(hashmaps),
        query_type_(query_type),
        start_level_(start_level),
        end_level_(end_level),
        inputs_(std::move(inputs)),
        start_index_(start_index),
        input_size_(input_size) {
    DCHECK_NE(input_size_, 0);
    DCHECK(!inputs_.empty());
  }

  inline const GraphBase* getDataGraph() const override { return graph_; }
  inline const std::vector<CandidateSetView>* getCandidates() const { return candidates_; }
  inline const CandidateHashmaps* getCandidateHashmaps() const { return candidate_hashmaps_; }

  inline const auto& getProfileInfo() const { return profile_info_; }
  inline const TaskStatus& getTaskStatus() const { return task_status_; }

  inline void setSuspendIntervalPtr(double* interval) { suspend_interval_ = interval; }
  inline void setSplitSize(uint32_t size) { split_size_ = size; }

  inline std::vector<std::pair<uint32_t, std::vector<CompressedSubgraphs>>> getSplits() { return std::move(splits_); }

  void run(uint32_t executor_idx) override;

  void profile(uint32_t executor_idx) override;

 protected:
  void logTask(uint64_t old_count, uint32_t executor_idx);

  /** @returns Whether to suspend task for creating new tasks on splits. */
  template <QueryType profile>
  bool splitInput();

  /** @returns True if success, otherwise task is pruned. */
  bool setupTraverseContexts();

  bool setupTraverseContextsWithoutCandidate();

  void setupOutputs();

  virtual std::unique_ptr<TraverseContext> createTraverseContext(const CandidateSetView* candidates,
                                                                 std::vector<CompressedSubgraphs>& outputs,
                                                                 uint32_t level, const TraverseOperator* op,
                                                                 QueryType profile,
                                                                 const unordered_set<VertexID>* candidate_hashmap) {
    return op->initTraverseContext(candidates, &outputs, graph_, profile, candidate_hashmap);
  }

  template <QueryType profile>
  bool execute(const std::vector<CompressedSubgraphs>& input, uint32_t start_index, uint32_t input_size, uint32_t level,
               uint32_t executor_idx);

  inline void checkSplits() {
    if (task_status_ == TaskStatus::Normal && !splits_.empty()) {
      task_status_ = TaskStatus::Split;
    }
  }
};

class TraverseTask : public TraverseChainTask {
 private:
  const ReorderedPartitionedGraph* partitioned_graph_;
  std::vector<CandidateScope> scopes_;
  std::vector<GraphView<GraphPartitionBase>> graph_views_;

 public:
  TraverseTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
               uint32_t batch_size, const std::vector<Operator*>* operators, std::unique_ptr<InputOperator>&& input_op,
               const std::vector<CandidateScope>& scopes, const GraphBase* graph,
               const std::vector<CandidateSetView>* candidates, QueryType query_type, uint32_t end_level,
               const CandidateHashmaps* hashmaps)
      : TraverseChainTask(query_id, task_id, stop_time, batch_size, operators, std::move(input_op), graph, candidates,
                          query_type, end_level, hashmaps),
        partitioned_graph_(dynamic_cast<const ReorderedPartitionedGraph*>(graph)),
        scopes_(scopes) {
    CHECK(partitioned_graph_ != nullptr);
    graph_views_ = setupGraphView(partitioned_graph_, *operators, scopes);
  }

  TraverseTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
               uint32_t batch_size, const std::vector<Operator*>* operators, const std::vector<CandidateScope>& scopes,
               const GraphBase* graph, const std::vector<CandidateSetView>* candidates, QueryType query_type,
               uint32_t start_level, uint32_t end_level, std::vector<CompressedSubgraphs>&& inputs,
               uint32_t start_index, uint32_t input_size, const CandidateHashmaps* hashmaps)
      : TraverseChainTask(query_id, task_id, stop_time, batch_size, operators, graph, candidates, query_type,
                          start_level, end_level, std::move(inputs), start_index, input_size, hashmaps),
        partitioned_graph_(dynamic_cast<const ReorderedPartitionedGraph*>(graph)),
        scopes_(scopes) {
    CHECK(partitioned_graph_ != nullptr);
    graph_views_ = setupGraphView(partitioned_graph_, *operators, scopes);
  }

  inline const std::vector<CandidateScope>& getScopes() const { return scopes_; }

 protected:
  std::unique_ptr<TraverseContext> createTraverseContext(const CandidateSetView* candidates,
                                                         std::vector<CompressedSubgraphs>& outputs, uint32_t level,
                                                         const TraverseOperator* op, QueryType profile,
                                                         const unordered_set<VertexID>* candidate_hashmap) override {
    return op->initTraverseContext(candidates, &outputs, &graph_views_[level - start_level_], profile,
                                   candidate_hashmap);
  }

  std::vector<GraphView<GraphPartitionBase>> setupGraphView(const ReorderedPartitionedGraph* g,
                                                            const std::vector<Operator*>& operators,
                                                            const std::vector<CandidateScope>& scopes);
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

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

// TODO(tatiana): move def to .cc
enum class TaskStatus : uint64_t { Normal, Suspended };
/**
 * For single-thread query execution
 */
class TraverseChainTask : public TaskBase {
 protected:
  const GraphBase* graph_ = nullptr;
  const uint32_t batch_size_;
  const std::vector<Operator*> operators_;
  const std::unique_ptr<InputOperator> input_op_ = nullptr;
  const std::vector<CandidateSetView>* candidates_ = nullptr;
  using CandidateHashmaps = std::vector<std::shared_ptr<unordered_set<VertexID>>>;
  const CandidateHashmaps* candidate_hashmaps_ = nullptr;

  std::vector<ProfileInfo> profile_info_;
  QueryType query_type_;
  std::vector<std::unique_ptr<TraverseContext>> traverse_context_;
  std::vector<std::vector<CompressedSubgraphs>> outputs_;
  uint32_t start_level_;
  uint32_t end_level_;
  TaskStatus task_status_;
  std::vector<CompressedSubgraphs> inputs_;
  uint32_t start_index_ = 0;
  uint32_t input_size_ = 0;

  std::chrono::time_point<std::chrono::steady_clock> start_time_;
  uint32_t suspended_level_ = 0;
  double suspend_interval_ = 0;
  uint32_t split_level_ = 0;
  uint32_t split_size_ = 1;
  std::vector<std::pair<uint32_t, std::vector<CompressedSubgraphs>>> splits_;

  VertexID seed_data_vertex_ = INVALID_VERTEX_ID;

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
        query_type_(query_type),
        start_level_(0),
        end_level_(operators_.size() - 1),
        task_status_(TaskStatus::Normal) {
    CHECK(!candidates_->empty());
  }

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
                    uint32_t batch_size, const std::vector<Operator*>& ops, std::unique_ptr<InputOperator>&& input_op,
                    const GraphBase* graph, const std::vector<CandidateSetView>* candidates, QueryType query_type,
                    uint32_t end_level, TaskStatus task_status, const CandidateHashmaps* hashmaps)
      : TaskBase(qid, tid, stop_time),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        input_op_(std::move(input_op)),
        candidates_(candidates),
        candidate_hashmaps_(hashmaps),
        query_type_(query_type),
        start_level_(0),
        end_level_(end_level),
        task_status_(task_status) {
    start_level_ = 0;
    if (end_level == 0) {
      end_level_ = ops.size() - 1;
    }
    CHECK(!candidates_->empty());
  }

  // for task starting from partial embeddings
  TraverseChainTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                    uint32_t batch_size, const std::vector<Operator*>& ops, const GraphBase* graph,
                    const std::vector<CandidateSetView>* candidates, QueryType query_type, uint32_t start_level,
                    uint32_t end_level, TaskStatus task_status, const std::vector<CompressedSubgraphs>&& inputs,
                    uint32_t start_index, uint32_t input_size, const CandidateHashmaps* hashmaps)
      : TaskBase(qid, tid, stop_time),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        candidates_(candidates),
        candidate_hashmaps_(hashmaps),
        query_type_(query_type),
        start_level_(start_level),
        end_level_(end_level),
        task_status_(task_status),
        inputs_(std::move(inputs)),
        start_index_(start_index),
        input_size_(input_size) {
    CHECK(!candidates_->empty());
    DCHECK_NE(input_size_, 0);
    DCHECK(!inputs_.empty());
  }

  const GraphBase* getDataGraph() const override { return graph_; }
  const std::vector<CandidateSetView>* getCandidates() const { return candidates_; }
  const auto& getProfileInfo() const { return profile_info_; }
  const TaskStatus& getTaskStatus() const { return task_status_; }
  const uint32_t getEndLevel() const { return end_level_; }
  std::vector<CompressedSubgraphs>& getLastOutput() { return outputs_[end_level_ - 1]; }
  const uint32_t getLastOutputSize() const { return traverse_context_[end_level_ - 1]->getOutputSize(); }

  void setSuspendInterval(double interval) { suspend_interval_ = interval; }
  void setSplitSize(uint32_t size) { split_size_ = size; }

  std::vector<std::pair<uint32_t, std::vector<CompressedSubgraphs>>> getSplits() { return std::move(splits_); }

  /** @returns Whether to suspend task for creating new tasks on splits. */
  template <QueryType profile>
  bool splitInput() {
    if (split_level_ <= start_level_) {  // if split from start level, directly get from input
      if (traverse_context_[0]->canSplitInput()) {
        auto[ptr, size] = traverse_context_[0]->splitInput();
        // LOG(INFO) << "Split start " << start_level_ << " size " << size << " interval " << suspend_interval_;
        splits_.emplace_back(start_level_, std::vector<CompressedSubgraphs>(ptr, ptr + size));
        if (splits_.size() == split_size_) {
          return true;
        }
      }
      split_level_ = start_level_ + 1;
    }
    auto split_max_level = std::min(suspended_level_, (uint32_t)operators_.size() - 2);
    if (split_level_ == split_max_level) return split_level_ >= operators_.size() - 2;

    if (traverse_context_[split_level_ - start_level_]->canSplitInput()) {
      auto[ptr, size] = traverse_context_[split_level_ - start_level_]->splitInput();
      // LOG(INFO) << "Split input " << split_level_ << " size " << size << " interval " << suspend_interval_;
      splits_.emplace_back(split_level_, std::vector<CompressedSubgraphs>(ptr, ptr + size));
      if (splits_.size() == split_size_) {
        return true;
      }
    }

    // split from intermediate level
    while (split_level_ < split_max_level && splits_.size() < split_size_) {
      auto ctx = traverse_context_[split_level_ - start_level_ - 1].get();
      auto original_buffer = ctx->getOutputs();
      std::vector<CompressedSubgraphs> temp_output_buffer = *original_buffer;
      auto original_size = ctx->getOutputSize();
      ctx->setOutputBuffer(temp_output_buffer, 0);
      auto traverse_op = dynamic_cast<TraverseOperator*>(operators_[split_level_ - 1]);
      uint32_t size = 0;
      if
        constexpr(isProfileMode(profile)) size = traverse_op->expandAndProfile(batch_size_, ctx);
      else
        size = traverse_op->expand(batch_size_, ctx);
      ctx->setOutputBuffer(*original_buffer, original_size);
      if (size == 0) {
        ++split_level_;
      } else {
        temp_output_buffer.erase(temp_output_buffer.begin() + size, temp_output_buffer.end());
        // LOG(INFO) << "Split middle " << split_level_ << " size " << size << " interval " << suspend_interval_;
        splits_.emplace_back(split_level_, std::move(temp_output_buffer));
      }
    }
    return splits_.size() == split_size_ || split_level_ >= operators_.size() - 2;
  }

  void run(uint32_t executor_idx) override {
    start_time_ = std::chrono::steady_clock::now();
    if (task_status_ == TaskStatus::Normal) {
      if (input_size_ == 0) {
        if (seed_data_vertex_ != INVALID_VERTEX_ID) {
          inputs_ = std::vector<CompressedSubgraphs>({CompressedSubgraphs(seed_data_vertex_)});
          input_size_ = 1;
        } else {
          if (verboseExecutionLog()) {
            std::stringstream ss;
            for (auto& x : *candidates_) {
              ss << x.size() << " ";
            }
            LOG(INFO) << "-------- Task " << task_id_ << " candidate sizes " << ss.str();
          }
          inputs_ = std::move(input_op_->getInputs(graph_, *candidates_));
          // if input is pruned to empty, return
          if (inputs_.empty()) {
            return;
          }
          input_size_ = inputs_.size();
        }
        LOG(INFO) << input_size_;
      }
      setupOutputs();
      if (seed_data_vertex_ != INVALID_VERTEX_ID) {
        if (!setupTraverseContextsWithoutCandidate()) return;
      } else {
        if (!setupTraverseContexts()) return;
      }
    }
    uint64_t old_count = 0;
    uint64_t new_count = 0;

    if (end_level_ == operators_.size() - 1) {
      old_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    }
    // TODO(engineering): support match limit
    // auto profile_output = "Task_" + std::to_string(task_id_);
    // ProfilerStart(profile_output.data());
    execute<QueryType::Execute>(inputs_, start_index_, input_size_, start_level_, executor_idx);
    // ProfilerStop();

    if (task_status_ == TaskStatus::Normal) {
      CHECK(splits_.empty());
    }
    if (end_level_ == operators_.size() - 1) {
      new_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    }
    auto end = std::chrono::steady_clock::now();
    if (shortExecutionLog()) {
      if (start_level_ == 0) {
        LOG(INFO) << "Task " << task_id_ << " input " << inputs_.size() << '/'
                  << getNumSubgraphs(inputs_, 0, inputs_.size()) << " count " << (new_count - old_count) << " ("
                  << old_count << " to " << new_count << "), time usage: " << toSeconds(start_time_, end) << "s.";
      } else {
        LOG(INFO) << "Inherit Task " << task_id_ << " input " << inputs_.size() << " count " << (new_count - old_count)
                  << " (" << old_count << " to " << new_count << "), time usage: " << toSeconds(start_time_, end)
                  << "s.";
      }
    }
  }

  void profile(uint32_t executor_idx) override {
    if (task_status_ == TaskStatus::Normal) {
      profile_info_.resize(operators_.size() + 1);  // traverse chain + input operator
      if (input_size_ == 0) {
        input_op_->inputAndProfile(graph_, *candidates_, &inputs_, &profile_info_.front());
        if (inputs_.empty()) return;
        input_size_ = inputs_.size();
      }
      setupOutputs();
      if (!setupTraverseContexts()) return;
    }
    // auto profile_output = "Profile_Task_" + std::to_string(task_id_);
    // ProfilerStart(profile_output.data());
    // TODO(engineering): support match limit
    start_time_ = std::chrono::steady_clock::now();
    uint64_t old_count = 0;
    uint64_t new_count = 0;
    if (end_level_ == operators_.size() - 1) {
      old_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    }
    execute<QueryType::Profile>(inputs_, start_index_, input_size_, start_level_, executor_idx);
    if (task_status_ == TaskStatus::Normal) {
      CHECK(splits_.empty());
      DCHECK_LT(start_level_, operators_.size());
      DCHECK_LT(end_level_, operators_.size());
      for (uint32_t i = start_level_; i < end_level_; ++i) {
        profile_info_[i + 1] += *traverse_context_[i - start_level_];
      }
    }
    if (end_level_ == operators_.size() - 1) {
      new_count = dynamic_cast<OutputOperator*>(operators_.back())->getOutput()->getCount(executor_idx);
    }
    auto end = std::chrono::steady_clock::now();
    if (shortExecutionLog()) {
      if (task_status_ == TaskStatus::Suspended) {
        LOG(INFO) << "Suspended Task " << task_id_ << " split " << split_level_ << " size " << splits_.size()
                  << " suspend " << suspended_level_ << '/' << (operators_.size() - 1)
                  << " time usage: " << toSeconds(start_time_, end) << "s. suspend interval " << suspend_interval_;
      } else if (start_level_ == 0) {
        LOG(INFO) << "Task " << task_id_ << " input " << inputs_.size() << '/'
                  << getNumSubgraphs(inputs_, 0, inputs_.size()) << " count " << (new_count - old_count) << " ("
                  << old_count << " to " << new_count << "), time usage: " << toSeconds(start_time_, end)
                  << "s. suspend interval " << suspend_interval_;
      } else {
        LOG(INFO) << "Inherit Task " << task_id_ << ':' << start_level_ << " input " << inputs_.size() << " count "
                  << (new_count - old_count) << " (" << old_count << " to " << new_count
                  << "), time usage: " << toSeconds(start_time_, end) << "s.";
      }
    }
    // ProfilerStop();
  }

 protected:
  /** @returns True if success, otherwise task is pruned. */
  bool setupTraverseContexts() {
    auto op = operators_[start_level_];
    auto traverse = dynamic_cast<TraverseOperator*>(op);
    uint32_t level = start_level_;
    traverse_context_.reserve(end_level_ - start_level_);
    while (level < end_level_) {
      DCHECK(traverse != nullptr);
      DCHECK_LT(traverse->getTargetQueryVertex(), candidates_->size());
      const unordered_set<VertexID>* hashmap = nullptr;
      if (candidate_hashmaps_ != nullptr) {
        DCHECK_LT(traverse->getTargetQueryVertex(), candidate_hashmaps_->size());
        hashmap = (*candidate_hashmaps_)[traverse->getTargetQueryVertex()].get();
      }
      // now assume all query vertices have candidate sets
      auto ctx = createTraverseContext(&(*candidates_)[traverse->getTargetQueryVertex()],
                                       outputs_[level - start_level_], level, traverse, query_type_, hashmap);
      if (ctx == nullptr) return false;
      traverse_context_.emplace_back(std::move(ctx));
      // LOG(INFO) << "set candidates for " << traverse->getTargetQueryVertex() << " " << traverse->toString();
      traverse = dynamic_cast<TraverseOperator*>(traverse->getNext());
      ++level;
    }
    return true;
  }

  bool setupTraverseContextsWithoutCandidate() {
    auto op = operators_[start_level_];
    auto traverse = dynamic_cast<TraverseOperator*>(op);
    uint32_t level = start_level_;
    traverse_context_.reserve(end_level_ - start_level_);
    while (level < end_level_) {
      DCHECK(traverse != nullptr);
      // now assume all query vertices have candidate sets
      auto ctx = createTraverseContext(nullptr, outputs_[level - start_level_], level, traverse, query_type_, nullptr);
      if (ctx == nullptr) return false;
      traverse_context_.emplace_back(std::move(ctx));
      // LOG(INFO) << "set candidates for " << traverse->getTargetQueryVertex() << " " << traverse->toString();
      traverse = dynamic_cast<TraverseOperator*>(traverse->getNext());
      ++level;
    }
    return true;
  }

  void setupOutputs() {
    uint32_t key_size = inputs_.front().getKeys().size();
    uint32_t set_size = inputs_.front().getSets().size();

    auto output_size = std::make_pair(key_size, key_size + set_size);
    outputs_.resize(end_level_ - start_level_);
    for (uint32_t i = start_level_; i < end_level_; ++i) {
      output_size = ((const TraverseOperator*)operators_[i])->getOutputSize(output_size);
      outputs_[i - start_level_].resize(batch_size_, CompressedSubgraphs(output_size.first, output_size.second));
    }
  }

  virtual std::unique_ptr<TraverseContext> createTraverseContext(const CandidateSetView* candidates,
                                                                 std::vector<CompressedSubgraphs>& outputs,
                                                                 uint32_t level, const TraverseOperator* op,
                                                                 QueryType profile,
                                                                 const unordered_set<VertexID>* candidate_hashmap) {
    return op->initTraverseContext(candidates, &outputs, graph_, profile, candidate_hashmap);
  }

  template <QueryType profile>
  bool execute(const std::vector<CompressedSubgraphs>& input, uint32_t start_index, uint32_t input_size, uint32_t level,
               uint32_t executor_idx) {
    if (isTimeOut()) {
      return true;
    }
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      if
        constexpr(isProfileMode(profile)) {
          return output_op->validateAndOutputAndProfile(input, 0, input_size, executor_idx, &profile_info_[level + 1]);
        }
      uint32_t start = 0;
      return output_op->validateAndOutput(input, start, input_size, executor_idx);
    }

    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    auto ctx = traverse_context_[level - start_level_].get();
    if (task_status_ == TaskStatus::Normal) {
      DCHECK_NOTNULL(ctx->getOutputs());
      ctx->setInput(input, start_index, start_index + input_size);
      if (suspend_interval_ != 0 && toSeconds(start_time_, std::chrono::steady_clock::now()) >= suspend_interval_) {
        suspended_level_ = level;
        if (splitInput<profile>()) {  // do not suspend if no splits
          task_status_ = TaskStatus::Suspended;
          return true;
        }
      }
    } else if (level == suspended_level_) {  // when a suspended task is resumed, we need to first recreate the stack
                                             // until the last level
      task_status_ = TaskStatus::Normal;
    }

    bool finished = false;
    while (true) {
      if (isTimeOut()) {
        return true;
      }
      uint32_t size = 0;
      if (task_status_ == TaskStatus::Normal) {
        ctx->resetOutputs();
        DCHECK_EQ(ctx->getOutputs()->size(), batch_size_) << "level " << level;
        if
          constexpr(isProfileMode(profile)) size = traverse_op->expandAndProfile(batch_size_, ctx);
        else
          size = traverse_op->expand(batch_size_, ctx);

        if (size == 0) {
          break;
        }
      }
      DCHECK_EQ(outputs_[level - start_level_].size(), batch_size_) << "level " << level;
      if (execute<profile>(outputs_[level - start_level_], 0, size, level + 1, executor_idx)) {
        finished = true;
        break;
      }
    }
    if
      constexpr(isProfileMode(profile)) {
        if (task_status_ == TaskStatus::Normal) {
          ctx->total_input_size += ctx->getTotalInputSize();
        }
      }

    return finished;
  }
};

class TraverseTask : public TraverseChainTask {
 private:
  const ReorderedPartitionedGraph* partitioned_graph_;
  std::vector<CandidateScope> scopes_;
  std::vector<GraphView<GraphPartitionBase>> graph_views_;

 public:
  TraverseTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
               uint32_t batch_size, const std::vector<Operator*>& operators, std::unique_ptr<InputOperator>&& input_op,
               const std::vector<CandidateScope>& scopes, const GraphBase* graph,
               const std::vector<CandidateSetView>* candidates, QueryType query_type, uint32_t end_level,
               TaskStatus task_status, const CandidateHashmaps* hashmaps)
      : TraverseChainTask(query_id, task_id, stop_time, batch_size, operators, std::move(input_op), graph, candidates,
                          query_type, end_level, task_status, hashmaps),
        partitioned_graph_(dynamic_cast<const ReorderedPartitionedGraph*>(graph)),
        scopes_(scopes) {
    CHECK(partitioned_graph_ != nullptr);
    graph_views_ = setupGraphView(partitioned_graph_, operators, scopes);
  }

  TraverseTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
               uint32_t batch_size, const std::vector<Operator*>& operators, const std::vector<CandidateScope>& scopes,
               const GraphBase* graph, const std::vector<CandidateSetView>* candidates, QueryType query_type,
               uint32_t start_level, uint32_t end_level, TaskStatus task_status,
               const std::vector<CompressedSubgraphs>&& inputs, uint32_t start_index, uint32_t input_size,
               const CandidateHashmaps* hashmaps)
      : TraverseChainTask(query_id, task_id, stop_time, batch_size, operators, graph, candidates, query_type,
                          start_level, end_level, task_status, std::move(inputs), start_index, input_size, hashmaps),
        partitioned_graph_(dynamic_cast<const ReorderedPartitionedGraph*>(graph)),
        scopes_(scopes) {
    CHECK(partitioned_graph_ != nullptr);
    graph_views_ = setupGraphView(partitioned_graph_, operators, scopes);
  }

  const CandidateHashmaps* getCandidateHashmaps() const { return candidate_hashmaps_; }

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
                                                            const std::vector<CandidateScope>& scopes) {
    std::vector<GraphView<GraphPartitionBase>> data_graphs_for_operators;
    CHECK_GT(operators.size(), 0);
    // size_t len = operators.size() - 1;
    size_t len = end_level_ - start_level_;
    data_graphs_for_operators.reserve(len);
    for (size_t i = start_level_; i < end_level_; ++i) {
      auto op = operators[i];
      auto traverse_op = dynamic_cast<TraverseOperator*>(op);
      CHECK(traverse_op != nullptr) << i << '/' << end_level_ << ' ' << op->toString();
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

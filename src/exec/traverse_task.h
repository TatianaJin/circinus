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

#include "exec/result.h"
#include "exec/task.h"
#include "graph/graph.h"
#include "graph/graph_partition.h"
#include "graph/graph_view.h"
#include "ops/input_operator.h"
#include "ops/traverse_operator.h"

namespace circinus {

class TraverseTask : public TaskBase {
 private:
  const uint32_t batch_size_;
  std::unique_ptr<InputOperator> input_op_;
  const ReorderedPartitionedGraph* graph_;
  std::vector<Operator*> operators_;
  std::vector<CandidateSetView> candidates_;
  std::vector<TraverseContext> traverse_context_;
  std::vector<GraphView<GraphPartition*>> graph_views_;
  QueryType query_type;

 public:
  TraverseTask(QueryId query_id, TaskId task_id, uint32_t batch_size, std::vector<Operator*>& operators,
               std::unique_ptr<InputOperator>&& input_op, const std::vector<CandidateScope>& scopes,
               const GraphBase* graph, CandidateResult* candidates)
      : TaskBase(query_id, task_id),
        input_op_(std::move(input_op)),
        operators_(operators),
        graph_(dynamic_cast<const ReorderedPartitionedGraph*>(graph)),
        batch_size_(batch_size) {
    CHECK(graph_ != nullptr);
    CHECK_EQ(candidates->getMergedCandidates()->size(), scopes.size());
    auto partitioned_candidates = dynamic_cast<PartitionedCandidateResult*>(candidates);
    CHECK(partitioned_candidates != nullptr);

    candidates_.reserve(scopes.size());
    for (uint32_t i = 0; i < scopes.size(); ++i) {
      candidates_.emplace_back(&partitioned_candidates->getMergedCandidates(i), scopes[i],
                               partitioned_candidates->getCandidatePartitionOffsets(i));
    }
    graph_views_ = setupGraphView(graph_, operators, scopes);
  }

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override {
    const auto& inputs = input_op_->getInputs(graph_, candidates_);
    uint32_t len = operators_.size() - 1;
    traverse_context_.reserve(len);
    for (uint32_t i = 0; i < len; ++i) {
      if (i == 0) {
        traverse_context_.emplace_back(TraverseContext(0, &inputs, &graph_views_[i]));
      } else {
        const auto& last_output = traverse_context_.back().getOutputs();
        traverse_context_.emplace_back(TraverseContext(0, last_output, &graph_views_[i]));
      }
    }
    execute();
  }

  void profile() override {
    const auto& inputs = input_op_->getInputs(graph_, candidates_);
    uint32_t len = operators_.size() - 1;
    traverse_context_.reserve(len);
    for (uint32_t i = 0; i < len; ++i) {
      if (i == 0) {
        traverse_context_.emplace_back(TraverseContext(0, &inputs, &graph_views_[i]));
      } else {
        const auto& last_output = traverse_context_.back().getOutputs();
        traverse_context_.emplace_back(TraverseContext(0, last_output, &graph_views_[i]));
      }
      traverse_context_.back().type = query_type;
    }
    profile_execute();
  }

 protected:
  bool profile_execute(uint32_t level = 0) {
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      return output_op->validateAndOutputAndProfile(*traverse_context_[level - 1].getOutputs(), 0);
    }
    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    traverse_context_[level].setInputEndIndex();
    while (true) {
      traverse_context_[level].clearOutputs();
      auto size = traverse_op->expandAndProfile(batch_size_, &traverse_context_[level]);
      if (size == 0) {
        break;
      }
      if (profile_execute(level + 1)) {
        return true;
      }
    }
    return false;
  }

  bool execute(uint32_t level = 0) {
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      return output_op->validateAndOutput(*traverse_context_[level - 1].getOutputs(), 0);
    }
    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    traverse_context_[level].setInputEndIndex();
    while (true) {
      traverse_context_[level].clearOutputs();
      auto size = traverse_op->expand(batch_size_, &traverse_context_[level]);
      if (size == 0) {
        break;
      }
      if (execute(level + 1)) {
        return true;
      }
    }
    return false;
  }

  static std::vector<GraphView<GraphPartition*>> setupGraphView(const ReorderedPartitionedGraph* g,
                                                                const std::vector<Operator*>& operators,
                                                                const std::vector<CandidateScope>& scopes) {
    std::vector<GraphView<GraphPartition*>> data_graphs_for_operators;
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
 * For single-thread query execution
 */
class TraverseChainTask : public TaskBase {
  const Graph* graph_;
  const uint32_t batch_size_;
  std::vector<Operator*> operators_;
  std::vector<CandidateSetView> candidates_;
  const uint32_t input_candidate_index_;
  const bool inputs_are_keys_;

 public:
  TraverseChainTask(QueryId qid, TaskId tid, uint32_t batch_size, std::vector<Operator*>& ops, const Graph* graph,
                    std::vector<CandidateSetView>&& candidates, uint32_t input_candidate_index, bool inputs_are_keys)
      : TaskBase(qid, tid),
        graph_(graph),
        batch_size_(batch_size),
        operators_(ops),
        candidates_(std::move(candidates)),
        input_candidate_index_(input_candidate_index),
        inputs_are_keys_(inputs_are_keys) {}

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override {
    auto op = operators_[0];
    auto traverse = dynamic_cast<TraverseOperator*>(op);
    while (traverse != nullptr) {
      DCHECK_LT(traverse->getTargetQueryVertex(), candidates_.size());
      // now assume all query vertices have candidate sets
      //
      traverse->setCandidateSets(&candidates_[traverse->getTargetQueryVertex()]);
      LOG(INFO) << "set candidates for " << traverse->getTargetQueryVertex() << " " << traverse->toString();
      op = op->getNext();
      traverse = dynamic_cast<TraverseOperator*>(op);
    }
    if (inputs_are_keys_) {
      execute(std::vector<CompressedSubgraphs>(candidates_[input_candidate_index_].begin(),
                                               candidates_[input_candidate_index_].end()));
    } else {
      std::vector<CompressedSubgraphs> input;
      input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(candidates_[input_candidate_index_])));
      execute(input);
    }
  }

 protected:
  bool execute(const std::vector<CompressedSubgraphs>& inputs, uint32_t level = 0) {
    std::vector<CompressedSubgraphs> outputs;
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      return output_op->validateAndOutput(inputs, 0);
    }
    auto traverse_op = dynamic_cast<TraverseOperator*>(op);

    while (true) {
      outputs.clear();
      auto size = traverse_op->expand(&outputs, FLAGS_batch_size);
      if (size == 0) {
        break;
      }
      if (execute(outputs, level + 1)) {
        return true;
      }
    }
    return false;
  }
};

/**
 * For multi-thread query execution
 */
class MatchingParallelTask : public TaskBase {
 public:
  MatchingParallelTask(QueryId qid, TaskId tid) : TaskBase(qid, tid) {}

  uint32_t getNextLevel() { return task_id_ + 1; }
  virtual TraverseOperator* getNextOperator();
  virtual std::vector<CompressedSubgraphs>& getOutputs();
};

class MatchingParallelInputTask : public MatchingParallelTask {
  const Graph* graph_;
  std::unique_ptr<InputOperator> input_op_;
  std::vector<CompressedSubgraphs> outputs_;
  std::vector<CandidateSetView> candidates_;

 public:
  MatchingParallelInputTask(QueryId qid, TaskId tid, const Graph* graph, std::unique_ptr<InputOperator>&& input_op,
                            std::vector<CandidateSetView>&& candidates)
      : MatchingParallelTask(qid, tid),
        graph_(graph),
        input_op_(std::move(input_op)),
        candidates_(std::move(candidates)) {}

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override { outputs_ = std::move(input_op_->getInputs(graph_, candidates_)); }

  std::vector<CompressedSubgraphs>& getOutputs() override { return outputs_; }
};

/**
 * For multi-thread query execution
 */
class MatchingParallelTraverseTask : public MatchingParallelTask {
  const Graph* graph_;
  TraverseOperator* traverse_op_;
  TraverseContext* traverse_ctx_;
  std::vector<CandidateSetView> candidates_;

 public:
  MatchingParallelTraverseTask(QueryId qid, TaskId tid, const Graph* graph, TraverseOperator* traverse_op,
                               std::vector<CandidateSetView>&& candidates, TraverseContext* traverse_ctx,
                               const std::vector<CompressedSubgraphs>& inputs, uint32_t input_index,
                               uint32_t input_end_index)
      : MatchingParallelTask(qid, tid),
        graph_(graph),
        traverse_op_(traverse_op),
        candidates_(std::move(candidates)),
        traverse_ctx_(traverse_ctx) {}

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override { traverse_op_.expand(UINT32_MAX, traverse_ctx_); }

  std::vector<CompressedSubgraphs>& getOutputs() override { return *(traverse_ctx_->outputs); }

  TraverseOperator* getNextOperator() override { return traverse_op_->getNext(); }
};

// TODO(byli) output task
class MatchingParallelOutputTask : public MatchingParallelTask {};

}  // namespace circinus

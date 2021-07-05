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
#include "ops/input_operator.h"
#include "ops/traverse_operator.h"
#include "plan/operator_tree.h"

namespace circinus {

class TraverseTask : public TaskBase {
 private:
  std::unique_ptr<InputOperator> input_op_;
  OperatorTree* op_tree_;  // not owned
  const ReorderedPartitionedGraph* graph_;
  const uint32_t batch_size_;
  std::vector<CandidateSetView> candidates_;
  std::vector<GraphView<GraphPartitionBase>> graph_views_;

 public:
  TraverseTask(QueryId query_id, TaskId task_id, uint32_t batch_size, OperatorTree& op_tree,
               std::unique_ptr<InputOperator>&& input_op, const std::vector<CandidateScope>& scopes,
               const GraphBase* graph, CandidateResult* candidates)
      : TaskBase(query_id, task_id),
        input_op_(std::move(input_op)),
        op_tree_(&op_tree),
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
    graph_views_ = setupGraphView(graph_, op_tree, scopes);
  }

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override {
    auto inputs = input_op_->getInputs(graph_, candidates_);
    // op_tree_->execute()
    // TODO(tatiana)
    // traverse_->input();
    // traverse_->expand();
  }

  void profile() override {
    // TODO(tatiana)
    // traverse_->inputAndProfile();
    // traverse_->expandAndProfile();
  }

 protected:
  static std::vector<GraphView<GraphPartitionBase>> setupGraphView(const ReorderedPartitionedGraph* g,
                                                                   const OperatorTree& op_tree,
                                                                   const std::vector<CandidateScope>& scopes) {
    std::vector<GraphView<GraphPartitionBase>> data_graphs_for_operators;
    size_t len = op_tree.getOperatorSize() - 1;
    data_graphs_for_operators.reserve(len);
    for (size_t i = 0; i < len; ++i) {
      auto op = op_tree.getOperator(i);
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
  OperatorTree* const op_tree_;  // not owned
  std::vector<CandidateSetView> candidates_;
  const uint32_t input_candidate_index_;
  const bool inputs_are_keys_;

 public:
  TraverseChainTask(QueryId qid, TaskId tid, uint32_t batch_size, OperatorTree& ops, const Graph* graph,
                    std::vector<CandidateSetView>&& candidates, uint32_t input_candidate_index, bool inputs_are_keys)
      : TaskBase(qid, tid),
        graph_(graph),
        batch_size_(batch_size),
        op_tree_(&ops),
        candidates_(std::move(candidates)),
        input_candidate_index_(input_candidate_index),
        inputs_are_keys_(inputs_are_keys) {}

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override {
    auto op = op_tree_->root();
    auto traverse = dynamic_cast<TraverseOperator*>(op);
    while (traverse != nullptr) {
      DCHECK_LT(traverse->getTargetQueryVertex(), candidates_.size());
      // now assume all query vertices have candidate sets
      traverse->setCandidateSets(&candidates_[traverse->getTargetQueryVertex()]);
      LOG(INFO) << "set candidates for " << traverse->getTargetQueryVertex() << " " << traverse->toString();
      op = op->getNext();
      traverse = dynamic_cast<TraverseOperator*>(op);
    }
    if (inputs_are_keys_) {
      op_tree_->execute(graph_, std::vector<CompressedSubgraphs>(candidates_[input_candidate_index_].begin(),
                                                                 candidates_[input_candidate_index_].end()));
    } else {
      std::vector<CompressedSubgraphs> input;
      input.emplace_back(std::make_shared<std::vector<VertexID>>(candidates_[input_candidate_index_].begin(),
                                                                 candidates_[input_candidate_index_].end()));
      op_tree_->execute(graph_, input);
    }
  }
};

/**
 * For multi-thread query execution
 */
class MatchingParallelTask : public TaskBase {
 public:
  MatchingParallelTask(QueryId qid, TaskId tid) : TaskBase(qid, tid) {}

  uint32_t getNextLevel() const { return getTaskId() + 1; }
  virtual TraverseOperator* getNextOperator() = 0;
  virtual std::vector<CompressedSubgraphs>& getOutputs() = 0;
};

class MatchingParallelInputTask : public MatchingParallelTask {
  const GraphBase* graph_ = nullptr;
  const InputOperator* input_op_ = nullptr;
  std::vector<CompressedSubgraphs> outputs_;
  std::vector<CandidateSetView> candidates_;

 public:
  MatchingParallelInputTask(QueryId qid, TaskId tid, const GraphBase* graph, const InputOperator* input_op,
                            std::vector<CandidateSetView>&& candidates)
      : MatchingParallelTask(qid, tid), graph_(graph), input_op_(input_op), candidates_(std::move(candidates)) {}

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override { outputs_ = input_op_->getInputs(graph_, candidates_); }

  std::vector<CompressedSubgraphs>& getOutputs() override { return outputs_; }

  TraverseOperator* getNextOperator() override { return dynamic_cast<TraverseOperator*>(input_op_->getNext()); }
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
        traverse_ctx_(traverse_ctx) {
    traverse_op_->input(inputs, input_index, input_end_index, graph_, traverse_ctx_);
  }

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override { traverse_op_->expand(UINT32_MAX, traverse_ctx_); }

  std::vector<CompressedSubgraphs>& getOutputs() override { return *(traverse_ctx_->outputs); }

  TraverseOperator* getNextOperator() override { return dynamic_cast<TraverseOperator*>(traverse_op_->getNext()); }
};

}  // namespace circinus

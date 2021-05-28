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

#include <vector>

#include "exec/task.h"
#include "graph/graph.h"
#include "ops/traverse_operator.h"
#include "plan/operator_tree.h"

namespace circinus {

class TraverseTask : public TaskBase {
 private:
  // TODO(tatiana): change to const functor?
  TraverseOperator* traverse_;  // not owned
  const Graph* graph_;

 public:
  TraverseTask(QueryId query_id, TaskId task_id, uint32_t shard_id, TraverseOperator* traverse, const Graph* graph)
      : TaskBase(query_id, task_id), traverse_(traverse), graph_(graph) {}

  void run() override {
    // TODO(tatiana)
    // traverse_->input();
    // traverse_->expand();
  }

  void profile() override {
    // TODO(tatiana)
    // traverse_->inputAndProfile();
    // traverse_->expandAndProfile();
  }
};

/**
 * For single-thread query execution
 */
class TraverseChainTask : public TaskBase {
  const Graph* graph_;
  const uint32_t batch_size_;
  OperatorTree* const op_tree_;                    // not owned
  std::vector<std::vector<VertexID>> candidates_;  // not owned
  const uint32_t input_candidate_index_;
  const bool inputs_are_keys_;

 public:
  TraverseChainTask(QueryId qid, TaskId tid, uint32_t batch_size, OperatorTree& ops, const Graph* graph,
                    std::vector<std::vector<std::vector<VertexID>>>& candidates, uint32_t input_candidate_index,
                    bool inputs_are_keys)
      : TaskBase(qid, tid),
        graph_(graph),
        batch_size_(batch_size),
        op_tree_(&ops),
        candidates_(candidates.size()),
        input_candidate_index_(input_candidate_index),
        inputs_are_keys_(inputs_are_keys) {
    for (uint32_t i = 0; i < candidates.size(); ++i) {
      CHECK_EQ(candidates[i].size(), 1);
      candidates_[i] = std::move(candidates[i].front());
    }
  }

  void run() override {
    auto op = op_tree_->root();
    auto traverse = dynamic_cast<TraverseOperator*>(op);
    while (traverse != nullptr) {
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
      input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(candidates_[input_candidate_index_])));
      op_tree_->execute(graph_, input);
    }
  }
};

}  // namespace circinus

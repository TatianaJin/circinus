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

#include "exec/task.h"
#include "graph/graph.h"
#include "ops/traverse_operator.h"

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

}  // namespace circinus

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

#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "utils/query_utils.h"

namespace circinus {

class TaskBase {
 private:
  uint16_t query_id_;
  uint16_t task_id_;

 public:
  TaskBase(QueryId query_id, TaskId task_id) : query_id_(query_id), task_id_(task_id) {}

  virtual ~TaskBase() {}

  bool isBefore(TaskBase* other) const {
    return (query_id_ == other->query_id_) ? task_id_ < other->task_id_ : query_id_ < other->query_id_;
  }

  inline auto getQueryId() const { return query_id_; }
  inline auto getTaskId() const { return task_id_; }

  virtual const GraphBase* getDataGraph() const = 0;
  virtual void run() = 0;
  virtual void profile() { run(); }
};

class Task {
  uint32_t level_;
  std::vector<CompressedSubgraphs> input_;
  const Graph* data_graph_;

 public:
  Task(uint32_t level, std::vector<CompressedSubgraphs>&& input, const Graph* graph)
      : level_(level), input_(std::move(input)), data_graph_(graph) {}

  inline uint32_t getLevel() const { return level_; }
  inline const auto& getInput() const { return input_; }
  inline const Graph* getDataGraph() const { return data_graph_; }
};

}  // namespace circinus

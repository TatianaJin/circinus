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

#include <chrono>
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class TaskBase {
 protected:
  uint16_t query_id_;
  uint16_t task_id_;
  double time_ = 0;
  std::chrono::time_point<std::chrono::steady_clock> stop_time_;

 public:
  TaskBase(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time)
      : query_id_(query_id), task_id_(task_id), stop_time_(stop_time) {}

  virtual ~TaskBase() {}

  bool isBefore(TaskBase* other) const {
    return (query_id_ == other->query_id_) ? task_id_ > other->task_id_ : query_id_ < other->query_id_;
  }

  inline bool isTimeOut() const { return std::chrono::steady_clock::now() >= stop_time_; }
  inline auto getQueryId() const { return query_id_; }
  inline auto getTaskId() const { return task_id_; }
  inline auto getStopTime() const { return stop_time_; }

  virtual const GraphBase* getDataGraph() const = 0;
  virtual void run(uint32_t executor_idx) = 0;
  virtual void profile(uint32_t executor_idx) { run(executor_idx); }

  inline void runWithTiming(uint32_t executor_idx) {
    auto start = std::chrono::high_resolution_clock::now();
    run(executor_idx);
    time_ += toSeconds(start, std::chrono::high_resolution_clock::now());
  }

  /* @returns Execution time in seconds */
  inline double getExecutionTime() const { return time_; }
};

}  // namespace circinus

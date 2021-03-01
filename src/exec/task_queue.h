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

#include <mutex>
#include <queue>
#include <vector>

#include "exec/task.h"
#include "graph/compressed_subgraphs.h"

namespace circinus {

/** Concurrent task priority queue prioritizing execution of deep tasks. */
class TaskQueue {
  struct TaskComparator {
    bool operator()(Task* a, Task* b) { return a->getLevel() < b->getLevel(); }
  };
  std::priority_queue<Task*, std::vector<Task*>, TaskComparator> task_queue_;
  std::mutex mu_;

 public:
  /** A task is newed */
  inline void putTask(uint32_t level, std::vector<CompressedSubgraphs>&& input, const Graph* graph) {
    std::lock_guard<std::mutex> lock(mu_);
    task_queue_.push(new Task(level, std::move(input), graph));
  }

  /** A task is newed */
  inline void putTaskUnSafe(uint32_t level, std::vector<CompressedSubgraphs>&& input, const Graph* graph) {
    task_queue_.push(new Task(level, std::move(input), graph));
  }

  /** The caller should delete the pointer after use */
  inline Task* getTask() {
    std::lock_guard<std::mutex> lock(mu_);
    if (task_queue_.empty()) return nullptr;
    auto task = task_queue_.top();
    task_queue_.pop();
    return task;
  }

  /** The caller should delete the pointer after use */
  inline Task* getTaskUnsafe() {
    if (task_queue_.empty()) return nullptr;
    auto task = task_queue_.top();
    task_queue_.pop();
    return task;
  }

  inline size_t size() {
    std::lock_guard<std::mutex> lock(mu_);
    return task_queue_.size();
  }
};

}  // namespace circinus

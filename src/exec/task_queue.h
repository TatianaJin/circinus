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

#include <condition_variable>
#include <mutex>
#include <queue>
#include <vector>

#include "exec/task.h"
#include "graph/compressed_subgraphs.h"

namespace circinus {

class ThreadsafeTaskQueue {
 private:
  struct Comparator {
    bool operator()(TaskBase* a, TaskBase* b) { return b->isBefore(a); }
  };

  // not storing unique_ptr because priority_queue does not support non-const top()
  std::priority_queue<TaskBase*, std::vector<TaskBase*>, Comparator> queue_;  // need to delete the pointers
  std::mutex mu_;
  std::condition_variable cv_;

 public:
  ~ThreadsafeTaskQueue() {
    while (!queue_.empty()) {
      delete queue_.top();
      queue_.pop();
    }
  }

  void putTask(TaskBase* task) {
    {
      std::lock_guard<std::mutex> lock(mu_);
      queue_.push(task);
    }
    cv_.notify_one();
  }

  inline void putTask(std::unique_ptr<TaskBase>&& task) { putTask(task.release()); }

  std::unique_ptr<TaskBase> getTask() {
    std::unique_lock<std::mutex> lock(mu_);
    cv_.wait(lock, [this]() { return !queue_.empty(); });
    auto ret = std::unique_ptr<TaskBase>(queue_.top());
    queue_.pop();
    return ret;
  }
};

/** Concurrent task priority queue prioritizing execution of deep tasks. */
class[[deprecated]] TaskQueue {
  struct TaskComparator {
    bool operator()(Task* a, Task* b) { return a->getLevel() < b->getLevel(); }
  };
  std::priority_queue<Task*, std::vector<Task*>, TaskComparator> task_queue_;
  std::mutex mu_;
  std::condition_variable cv_;
  uint32_t n_subscribers_ = 0;
  uint32_t wait_count_ = 0;

 public:
  explicit TaskQueue(uint32_t n_subscribers) : n_subscribers_(n_subscribers) {}

  /** A task is newed */
  inline void putTask(uint32_t level, std::vector<CompressedSubgraphs> && input, const Graph* graph) {
    bool notify = false;
    {
      std::lock_guard<std::mutex> lock(mu_);
      task_queue_.push(new Task(level, std::move(input), graph));
      notify = (task_queue_.size() == 1);
    }
    if (notify) {
      cv_.notify_one();
    }
  }

  /** A task is newed */
  inline void putTaskUnSafe(uint32_t level, std::vector<CompressedSubgraphs> && input, const Graph* graph) {
    task_queue_.push(new Task(level, std::move(input), graph));
  }

  /** The caller should delete the pointer after use */
  inline Task* getTask(uint32_t subscriber_id) {
    std::unique_lock<std::mutex> lock(mu_);
    // std::lock_guard<std::mutex> lock(mu_);
    if (task_queue_.empty()) {
      if (++wait_count_ == n_subscribers_) {  // all tasks finished
        LOG(INFO) << "tasks all done " << wait_count_ << "/" << n_subscribers_;
        cv_.notify_all();
        return nullptr;
      }
      while (task_queue_.empty()) {
        cv_.wait(lock);
        if (wait_count_ == n_subscribers_) return nullptr;
      }
      --wait_count_;
    }
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

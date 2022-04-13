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
#include <memory>
#include <mutex>
#include <queue>
#include <utility>
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
  bool shut_down_ = false;

 public:
  ~ThreadsafeTaskQueue() {
    while (!queue_.empty()) {
      delete queue_.top();
      queue_.pop();
    }
  }

  std::uint64_t getSize() { return queue_.size(); }

  void shutDown() {
    {
      std::lock_guard<std::mutex> lock(mu_);
      shut_down_ = true;
    }
    cv_.notify_all();
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
    cv_.wait(lock, [this]() { return shut_down_ || !queue_.empty(); });
    if (shut_down_) return nullptr;
    auto ret = std::unique_ptr<TaskBase>(queue_.top());
    queue_.pop();
    return ret;
  }
};

}  // namespace circinus

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

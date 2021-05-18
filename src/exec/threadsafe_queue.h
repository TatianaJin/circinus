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

namespace circinus {

/** FIFO multi thread queue */
template <typename T>
class ThreadsafeQueue {
 public:
  ThreadsafeQueue() = default;
  ~ThreadsafeQueue() = default;
  ThreadsafeQueue(const ThreadsafeQueue&) = delete;
  ThreadsafeQueue& operator=(const ThreadsafeQueue&) = delete;
  ThreadsafeQueue(ThreadsafeQueue&&) = delete;
  ThreadsafeQueue& operator=(ThreadsafeQueue&&) = delete;

  void Push(const T& elem_ptr) {
    mu_.lock();
    queue_.push(elem_ptr);
    mu_.unlock();
    cond_.notify_all();
  }

  void Push(T&& elem_ptr) {
    mu_.lock();
    queue_.push(std::move(elem_ptr));
    mu_.unlock();
    cond_.notify_all();
  }

  T WaitAndPop() {
    std::unique_lock<std::mutex> lk(mu_);
    cond_.wait(lk, [this] { return !queue_.empty(); });
    T elem = std::move(queue_.front());
    queue_.pop();
    return elem;
  }

 private:
  std::mutex mu_;
  std::queue<T> queue_;
  std::condition_variable cond_;
};

}  // namespace circinus

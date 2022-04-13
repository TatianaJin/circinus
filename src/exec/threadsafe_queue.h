#pragma once

#include <condition_variable>
#include <mutex>
#include <queue>
#include <utility>

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

  void push(const T& elem_ptr) {
    mu_.lock();
    queue_.push(elem_ptr);
    mu_.unlock();
    cond_.notify_all();
  }

  void push(T&& elem_ptr) {
    mu_.lock();
    queue_.push(std::move(elem_ptr));
    mu_.unlock();
    cond_.notify_all();
  }

  T waitAndPop() {
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

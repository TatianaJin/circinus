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

#include <thread>
#include <vector>

#include "exec/task_queue.h"
#include "graph/compressed_subgraphs.h"
#include "plan/execution_plan.h"

namespace circinus {

class ThreadPool {
  const ExecutionPlan* plan_;
  TaskQueue task_queue_;
  std::vector<std::thread> pool_;
  uint32_t n_threads_;
  std::vector<uint64_t> output_count_;

 public:
  ThreadPool(uint32_t n_threads, const ExecutionPlan* plan) : plan_(plan), task_queue_(), n_threads_(n_threads) {}

  ~ThreadPool() {
    for (auto& t : pool_) {
      t.join();
    }
    uint64_t count = 0;
    for (auto c : output_count_) {
      count += c;
    }
    LOG(INFO) << "find matches " << count;
  }

  inline void addInitTask(uint32_t level, std::vector<CompressedSubgraphs>&& input, const Graph* data_graph) {
    task_queue_.putTaskUnSafe(level, std::move(input), data_graph);
  }

  void start() {
    pool_.reserve(n_threads_);
    output_count_.resize(n_threads_, 0);
    pool_.push_back(std::thread([this]() {
      auto& operators = plan_->getOperators();
      while (true) {
        auto task = task_queue_.getTask();
        if (task == nullptr) break;
        operators.handleTask(task, &task_queue_, &output_count_[0]);
        delete task;
      }
    }));
    for (uint32_t i = 1; i < n_threads_; ++i) {
      pool_.push_back(std::thread([this, i]() {
        auto operators = plan_->cloneOperators();
        while (true) {
          auto task = task_queue_.getTask();
          if (task == nullptr) break;
          operators.handleTask(task, &task_queue_, &output_count_[i]);
          delete task;
        }
      }));
    }
  }
};

}  // namespace circinus

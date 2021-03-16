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
  ExecutionPlan* plan_;
  TaskQueue task_queue_;
  std::vector<std::thread> pool_;
  uint32_t n_threads_;

 public:
  ThreadPool(uint32_t n_threads, ExecutionPlan* plan) : plan_(plan), task_queue_(n_threads), n_threads_(n_threads) {}

  ~ThreadPool() {
    for (auto& t : pool_) {
      t.join();
    }
  }

  inline void addInitTask(uint32_t level, std::vector<CompressedSubgraphs>&& input, const Graph* data_graph) {
    task_queue_.putTaskUnSafe(level, std::move(input), data_graph);
  }

  void start() {
    pool_.reserve(n_threads_);
    pool_.push_back(std::thread([this]() {
      auto& operators = plan_->getOperators();
      bool stop = false;
      while (!stop) {
        auto task = task_queue_.getTask(0);
        if (task == nullptr) break;
        // LOG(INFO) << "input no. CompressedSubgraphs " << task->getInput().size() << ' ' << task->getLevel() << ' '
        // << operators.getOperator(task->getLevel())->toString();
        stop = operators.handleTask(task, &task_queue_, 0);
        delete task;
      }
    }));
    for (uint32_t i = 1; i < n_threads_; ++i) {
      pool_.push_back(std::thread([this, i]() {
        auto operators = plan_->cloneOperators();
        bool stop = false;
        while (!stop) {
          auto task = task_queue_.getTask(i);
          if (task == nullptr) {
            break;
          }
          stop = operators.handleTask(task, &task_queue_, i);
          delete task;
        }
      }));
    }
  }
};

}  // namespace circinus

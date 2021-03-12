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

#include <vector>

#include "ops/operator.h"

namespace circinus {

class Task;
class TaskQueue;

/** Uses flag: FLAGS_batch_size */
class OperatorTree {
  std::vector<Operator*> operators_;

 public:
  ~OperatorTree() { clear(); }

  void clear() {
    for (auto op : operators_) {
      delete op;
    }
    operators_.clear();
  }
  inline bool empty() const { return operators_.empty(); }

  inline Operator* root() const { return operators_.front(); }
  inline void push_back(Operator* t) { operators_.push_back(t); }
  inline void reserve(uint32_t size) { operators_.reserve(size); }

  inline OperatorTree clone() const {
    OperatorTree ret;
    ret.operators_.reserve(operators_.size());
    for (auto op : operators_) {
      ret.operators_.push_back(op->clone());
    }
    return ret;
  }

  bool handleTask(Task* task, TaskQueue* queue, uint32_t thread_id);
};

}  // namespace circinus

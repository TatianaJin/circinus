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

#include <utility>

#include "exec/task.h"

namespace circinus {

class ProfileTaskBase {
 public:
  // FIXME(tatiana)
  virtual ProfileInfo* getProfileInfo() const { return nullptr; }
};

template <typename T>
class ProfileTask : public T, public ProfileTaskBase {
  static_assert(std::is_base_of<TaskBase, T>::value, "T should inherit from TaskBase");

 public:
  template <typename... Args>
  ProfileTask(QueryId qid, TaskId tid, Args... args) : T(qid, tid, std::move(args)...) {}

  void run(uint32_t executor_idx) override { T::profile(executor_idx); }

  const GraphBase* getDataGraph() const override { return T::getDataGraph(); }
};

}  // namespace circinus

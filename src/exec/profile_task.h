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

#include "exec/task.h"

namespace circinus {

class ProfileTaskBase {
 public:
  // FIXME(tatiana)
};

template <typename T>
class ProfileTask : public TaskBase, public ProfileTaskBase {
  static_assert(std::is_base_of<TaskBase, T>::value, "T should inherit from TaskBase");
  T task_;

 public:
  template <typename... Args>
  ProfileTask(Args... args) : task_(args...) {}

  void run() override { task_.profile(); }

  const GraphBase* getDataGraph() const override { return task_.getDataGraph(); }
};

}  // namespace circinus

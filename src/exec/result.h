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

#include <memory>
#include <vector>

#include "exec/task.h"
#include "graph/types.h"
#include "utils/query_utils.h"

namespace circinus {

class Result {
 public:
  static std::unique_ptr<Result> newCandidateResult(TaskId n_tasks);
  static std::unique_ptr<Result> newExecutionResult();

  virtual void collect(TaskBase* task) = 0;
};

class CandidateResult : public Result {
 private:
  std::vector<std::vector<std::vector<VertexID>>> candidates_;  // {query vertex, {shard, candidates}}

 public:
  explicit CandidateResult(TaskId n_tasks) : candidates_(n_tasks) {}

  void collect(TaskBase* task) override;

  const std::vector<std::vector<std::vector<VertexID>>>& getCandidates() { return candidates_; }

  std::vector<VertexID> getCandidateCardinality() const {
    std::vector<VertexID> ret(candidates_.size(), 0);
    for (uint32_t i = 0; i < candidates_.size(); ++i) {
      for (auto& shard : candidates_[i]) {
        ret[i] += shard.size();
      }
    }
    return ret;
  }
};

class ExecutionResult : public Result {
 public:
  void collect(TaskBase* task) override;
  void* data() {
    // TODO(tatiana)
    return nullptr;
  }
};

}  // namespace circinus

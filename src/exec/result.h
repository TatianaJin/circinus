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
#include "utils/utils.h"

namespace circinus {

class Result {
 public:
  static std::unique_ptr<Result> newCandidateResult(TaskId n_tasks);
  static std::unique_ptr<Result> newExecutionResult();

  virtual ~Result() {}

  virtual void collect(TaskBase* task) = 0;
};

class CandidateResult : public Result {
 private:
  std::vector<std::vector<std::vector<VertexID>>> candidates_;  // {query vertex, {shard, candidates}}
  std::vector<std::vector<VertexID>> merged_candidates_;

 public:
  explicit CandidateResult(TaskId n_tasks) : candidates_(n_tasks) {}

  void collect(TaskBase* task) override;

  void merge();

  void remove_invalid(QueryVertexID query_vertex);

  std::vector<std::vector<VertexID>>* getMergedCandidates() { return &merged_candidates_; }

  const std::vector<std::vector<VertexID>>& getCandidates() const { return merged_candidates_; }
  std::vector<std::vector<VertexID>>& getCandidates() { return merged_candidates_; }

  std::vector<VertexID> getMergedCandidateCardinality() const {
    std::vector<VertexID> ret(merged_candidates_.size(), 0);
    for (QueryVertexID u = 0; u < merged_candidates_.size(); ++u) {
      LOG(INFO) << "QueryVertexID " << u << " candidates size " << merged_candidates_[u].size();
      ret[u] = merged_candidates_[u].size();
    }
    return ret;
  }

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
  uint64_t count_;

 public:
  void collect(TaskBase* task) override;

  void setCount(uint64_t count) { count_ = count; }

  // TODO(tatiana): now we only consider count as output
  void* data() { return &count_; }
};

}  // namespace circinus

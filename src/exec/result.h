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
  static std::unique_ptr<Result> newPartitionedCandidateResult(TaskId n_tasks, uint32_t n_partitions);
  static std::unique_ptr<Result> newExecutionResult();

  virtual ~Result() {}

  // FIXME(tatiana): not a suitable common interface?
  virtual void collect(TaskBase* task) = 0;
};

class CandidateResult : public Result {
 protected:
  std::vector<std::vector<std::vector<VertexID>>> candidates_;  // {query vertex, {shard, candidates}}
  std::vector<std::vector<VertexID>> merged_candidates_;

 public:
  explicit CandidateResult(TaskId n_tasks) : candidates_(n_tasks) {}

  virtual ~CandidateResult() {}

  void collect(TaskBase* task) override;

  virtual std::vector<std::vector<VertexID>> getCandidateCardinality() const;

  void merge();

  void remove_invalid(QueryVertexID query_vertex);

  std::vector<std::vector<VertexID>>* getMergedCandidates() { return &merged_candidates_; }
};

class PartitionedCandidateResult : public CandidateResult {
 public:
  explicit PartitionedCandidateResult(uint32_t n_qvs, uint32_t n_partitions) : CandidateResult(n_qvs) {
    for (auto& shards : candidates_) {
      shards.resize(n_partitions);
    }
  }

  void collect(TaskBase* task) override;

  std::vector<std::vector<VertexID>> getCandidateCardinality() const override;
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

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
#include "ops/output_operator.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class Result {
 public:
  static std::unique_ptr<Result> newCandidateResult(TaskId n_tasks);
  static std::unique_ptr<Result> newPartitionedCandidateResult(TaskId n_tasks, uint32_t n_partitions);
  static std::unique_ptr<Result> newExecutionResult();

  virtual ~Result() {}
};

class CandidateResult : public Result {
 protected:
  std::vector<std::vector<std::vector<VertexID>>> candidates_;  // {query vertex, {shard, candidates}}
  std::vector<std::vector<VertexID>> merged_candidates_;

 public:
  explicit CandidateResult(TaskId n_tasks) : candidates_(n_tasks), merged_candidates_(n_tasks) {}

  virtual ~CandidateResult() {}

  virtual void collect(TaskBase* task);

  virtual std::vector<std::vector<VertexID>> getCandidateCardinality() const;

  virtual void merge(TaskBase* task);

  virtual void removeInvalid(QueryVertexID query_vertex);

  std::vector<std::vector<VertexID>>* getMergedCandidates() { return &merged_candidates_; }

  const std::vector<VertexID>& getMergedCandidates(uint32_t idx) const { return merged_candidates_[idx]; }
};

class PartitionedCandidateResult : public CandidateResult {
  std::vector<std::vector<VertexID>> candidate_partition_offsets_;
  std::vector<std::vector<VertexID>> per_partition_candidate_cardinality_;

 public:
  explicit PartitionedCandidateResult(uint32_t n_qvs, uint32_t n_partitions)
      : candidate_partition_offsets_(n_qvs),
        per_partition_candidate_cardinality_(n_partitions, std::vector<VertexID>(n_qvs)),
        CandidateResult(n_qvs) {
    for (auto& shards : candidates_) {
      shards.resize(n_partitions);
    }
  }

  void collect(TaskBase* task) override;

  void merge(TaskBase* task) override;

  void removeInvalid(QueryVertexID query_vertex) override;

  std::vector<std::vector<VertexID>> getCandidateCardinality() const override {
    return per_partition_candidate_cardinality_;
  }

  inline const auto& getCandidatePartitionOffsets(uint32_t idx) const {
    CHECK(!candidate_partition_offsets_.empty()) << "merge() must be called before getCandidatePartitionOffsets";
    return candidate_partition_offsets_[idx];
  }
};

class ExecutionResult : public Result {
  uint64_t count_ = 0;  // FIXME(tatiana): remove
  std::vector<std::vector<CompressedSubgraphs>> inputs_;
  Outputs outputs_;

 public:
  // FIXME(tatiana): remove
  void setCount(uint64_t count) { count_ = count; }

  Outputs& getOutputs() { return outputs_; }

  // TODO(tatiana): now we only consider count as output
  void* data() { return &count_; }
};

}  // namespace circinus

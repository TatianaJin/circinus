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

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "exec/task.h"
#include "graph/candidate_set_view.h"
#include "graph/types.h"
#include "ops/output_operator.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class Result {
 public:
  static std::unique_ptr<Result> newCandidateResult(TaskId n_tasks);
  static std::unique_ptr<Result> newPartitionedCandidateResult(TaskId n_tasks, uint32_t n_partitions);
  static std::unique_ptr<Result> newExecutionResult(bool profile = false);

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

  inline std::vector<CandidateSetView> getCandidates() const {
    std::vector<CandidateSetView> res(candidates_.size());
    for (uint32_t i = 0; i < res.size(); ++i) {
      if (!merged_candidates_[i].empty()) {
        res[i].addRange(merged_candidates_[i].data(), merged_candidates_[i].data() + merged_candidates_[i].size());
      } else {
        std::vector<uint32_t> order(candidates_[i].size());
        std::iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(),
             [&](uint32_t l, uint32_t r) { return candidates_[i][l].front() <= candidates_[i][r].front(); });
        for (auto j : order) {
          auto& shard = candidates_[i][j];
          res[i].addRange(shard.data(), shard.data() + shard.size());
        }
      }
    }
    return res;
  }

  std::vector<std::vector<VertexID>>* getMergedCandidates() { return &merged_candidates_; }

  const std::vector<VertexID>& getMergedCandidates(uint32_t idx) const { return merged_candidates_[idx]; }
};

class PartitionedCandidateResult : public CandidateResult {
  std::vector<std::vector<VertexID>> candidate_partition_offsets_;
  std::vector<std::vector<VertexID>> per_partition_candidate_cardinality_;

 public:
  explicit PartitionedCandidateResult(uint32_t n_qvs, uint32_t n_partitions)
      : CandidateResult(n_qvs),
        candidate_partition_offsets_(n_qvs),
        per_partition_candidate_cardinality_(n_partitions, std::vector<VertexID>(n_qvs)) {
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
  QueryResult result_;
  std::vector<std::vector<CompressedSubgraphs>> inputs_;
  Outputs outputs_;

 public:
  virtual ~ExecutionResult() {}

  inline void setCount() { result_.embedding_count = outputs_.getCount(); }
  inline void addEnumerateTime(double time) { result_.enumerate_time += time; }
  inline void setElapsedExecutionTime(double time) { result_.elapsed_execution_time = time; }
  inline void setMatchingOrder(std::string&& order) { result_.matching_order = std::move(order); }

  inline Outputs& getOutputs() { return outputs_; }
  inline QueryResult& getQueryResult() { return result_; }

  virtual void collect(TaskBase* task) {}
};

class ProfiledExecutionResult : public ExecutionResult {
  std::vector<ProfileInfo> profiles_;
  std::vector<std::string> profiled_plan_str_;

 public:
  ProfiledExecutionResult(uint32_t n_profiles) : profiles_(n_profiles) {}

  void setProfiledPlan(const std::vector<Operator*>& ops) {
    CHECK_EQ(ops.size(), profiles_.size());
    auto size = ops.size();
    profiled_plan_str_.resize(size);
    for (uint32_t i = 0; i < size; ++i) {
      // FIXME(tatiana): populate profiled_plan_str_
    }
  }

  const auto& getProfiledPlanStrings() const { return profiled_plan_str_; }

  void collect(TaskBase* task) override {
    // FIXME(tatiana)
    // dynamic_cast<ProfileTaskBase*>(task)
  }
};

}  // namespace circinus

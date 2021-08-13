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

#include "exec/result.h"

#include <algorithm>
#include <memory>
#include <utility>

#include "exec/scan_task.h"
#include "exec/task.h"
#include "exec/traverse_task.h"
#include "plan/execution_plan.h"
#include "utils/utils.h"

namespace circinus {

std::unique_ptr<Result> Result::newCandidateResult(TaskId n_tasks) {
  return std::make_unique<CandidateResult>(n_tasks);
}

std::unique_ptr<Result> Result::newPartitionedCandidateResult(TaskId n_tasks, uint32_t n_partitions) {
  return std::make_unique<PartitionedCandidateResult>(n_tasks, n_partitions);
}

std::unique_ptr<Result> Result::newExecutionResult(bool profile, const std::vector<ExecutionPlan*>& plans) {
  if (profile) {
    std::vector<uint32_t> plan_sizes;
    plan_sizes.reserve(plans.size());
    for (auto& plan : plans) {
      plan_sizes.push_back(plan->getOperators().size() + 1);
    }
    return std::make_unique<ProfiledExecutionResult>(plan_sizes);
  }
  return std::make_unique<ExecutionResult>();
}

void CandidateResult::collect(TaskBase* task) {
  auto scan = dynamic_cast<ScanTask*>(task);
  auto& shard_candidates = scan->getScanContext().candidates;
  if (!shard_candidates.empty()) {
    candidates_[task->getTaskId()].push_back(std::move(shard_candidates));
  }
}

void CandidateResult::merge(TaskBase* task) {
  auto task_id = task->getTaskId();
  std::vector<uint32_t> order(candidates_[task_id].size());
  std::iota(order.begin(), order.end(), 0);
  sort(order.begin(), order.end(), [&](uint32_t l, uint32_t r) {
    if (candidates_[task_id][l].empty() || candidates_[task_id][r].empty()) {
      return true;
    }
    return candidates_[task_id][l].front() < candidates_[task_id][r].front();
  });

  for (uint32_t j : order) {
    merged_candidates_[task_id].insert(merged_candidates_[task_id].end(), candidates_[task_id][j].begin(),
                                       candidates_[task_id][j].end());
    candidates_[task_id][j].clear();
  }
}

void CandidateResult::removeInvalid(QueryVertexID query_vertex) {
  merged_candidates_[query_vertex].erase(
      std::remove_if(merged_candidates_[query_vertex].begin(), merged_candidates_[query_vertex].end(),
                     [invalid = INVALID_VERTEX_ID](VertexID vid) { return vid == invalid; }),
      merged_candidates_[query_vertex].end());
}

std::vector<std::vector<VertexID>> CandidateResult::getCandidateCardinality() const {
  std::vector<VertexID> ret(candidates_.size(), 0);
  for (uint32_t i = 0; i < candidates_.size(); ++i) {
    if (!merged_candidates_[i].empty()) {
      ret[i] = merged_candidates_[i].size();
    } else {
      for (auto& shard : candidates_[i]) {
        ret[i] += shard.size();
      }
    }
  }
  return {std::move(ret)};
}

void PartitionedCandidateResult::collect(TaskBase* task) {
  auto scan = dynamic_cast<ScanTask*>(task);
  uint32_t task_id = scan->getTaskId();
  uint32_t partition = scan->getPartition();
  candidates_[task_id][partition] = std::move(scan->getScanContext().candidates);
  per_partition_candidate_cardinality_[partition][task_id] = candidates_[task_id][partition].size();
}

void PartitionedCandidateResult::merge(TaskBase* task) {
  uint32_t task_id = task->getTaskId();
  candidate_partition_offsets_[task_id].resize(candidates_[task_id].size() + 1, 0);
  for (uint32_t j = 0; j < candidates_[task_id].size(); ++j) {
    candidate_partition_offsets_[task_id][j + 1] =
        candidate_partition_offsets_[task_id][j] + candidates_[task_id][j].size();
  }

  for (uint32_t j = 0; j < candidates_[task_id].size(); ++j) {
    merged_candidates_[task_id].insert(merged_candidates_[task_id].end(), candidates_[task_id][j].begin(),
                                       candidates_[task_id][j].end());
    candidates_[task_id][j].clear();
  }
}

void PartitionedCandidateResult::removeInvalid(QueryVertexID query_vertex) {
  VertexID last_invalid_sum = 0;
  VertexID invalid_sum = 0;
  VertexID valid_idx = 0;
  VertexID idx = 0;
  for (uint32_t i = 1; i < candidate_partition_offsets_[query_vertex].size(); ++i) {
    VertexID& offset = candidate_partition_offsets_[query_vertex][i];
    for (; idx < offset; ++idx) {
      if (merged_candidates_[query_vertex][idx] == INVALID_VERTEX_ID) {
        invalid_sum++;
      } else if (valid_idx++ != idx) {
        merged_candidates_[query_vertex][valid_idx - 1] = merged_candidates_[query_vertex][idx];
      }
    }
    CHECK_LE(invalid_sum, offset) << "Error: invalid vertex number greater than offset.";
    CHECK_LE(invalid_sum - last_invalid_sum, per_partition_candidate_cardinality_[i - 1][query_vertex])
        << "Error: per partition invalid vertex number greater than cardinality.";

    offset -= invalid_sum;
    per_partition_candidate_cardinality_[i - 1][query_vertex] -= invalid_sum - last_invalid_sum;
    last_invalid_sum = invalid_sum;
  }
  merged_candidates_[query_vertex].resize(valid_idx);
}

void ProfiledExecutionResult::collect(TaskBase* task) {
  auto traverse_task = dynamic_cast<TraverseChainTask*>(task);
  if (traverse_task != nullptr) {
    DCHECK_LT(traverse_task->getTaskId(), profiles_.size());
    auto& profile = profiles_[traverse_task->getTaskId()];
    auto size = profile.size();
    DCHECK_EQ(size, traverse_task->getProfileInfo().size());
    for (uint32_t i = 0; i < size; ++i) {
      profile[i] += traverse_task->getProfileInfo()[i];
    }
  } else {
    auto matching_parallel_task = dynamic_cast<MatchingParallelTask*>(task);
    CHECK(matching_parallel_task != nullptr);
    // TODO(profile): record intersection count in input op? now skip input task
    matching_parallel_task->collectProfileInfo(profiles_.front()[matching_parallel_task->getTaskId()]);
  }
}

void ProfiledExecutionResult::setProfiledPlan(uint32_t profile_idx, const std::vector<Operator*>& ops,
                                              const InputOperator* input_op) {
  CHECK_EQ(1 + ops.size(), profiles_[profile_idx].size());
  auto size = ops.size();
  profiled_plan_str_[profile_idx].reserve(size + 1);
  profiled_plan_str_[profile_idx].push_back(input_op->toProfileString(profiles_[profile_idx].front()));
  for (uint32_t i = 0; i < size; ++i) {
    profiled_plan_str_[profile_idx].push_back(ops[i]->toProfileString(profiles_[profile_idx][i + 1]));
  }
}

}  // namespace circinus

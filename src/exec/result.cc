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
#include "utils/utils.h"

namespace circinus {

std::unique_ptr<Result> Result::newCandidateResult(TaskId n_tasks) {
  return std::make_unique<CandidateResult>(n_tasks);
}

std::unique_ptr<Result> Result::newPartitionedCandidateResult(TaskId n_tasks, uint32_t n_partitions) {
  return std::make_unique<PartitionedCandidateResult>(n_tasks, n_partitions);
}

std::unique_ptr<Result> Result::newExecutionResult() { return std::make_unique<ExecutionResult>(); }

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
    for (auto& shard : candidates_[i]) {
      ret[i] += shard.size();
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
  CandidateResult::merge(task);
}

void PartitionedCandidateResult::removeInvalid(QueryVertexID query_vertex) {
  VertexID last_invalid_sum = 0;
  VertexID invalid_sum = 0;
  VertexID valid_idx = 0;
  for (uint32_t i = 1, j = 0; i < candidate_partition_offsets_[query_vertex].size(); ++i) {
    VertexID& offset = candidate_partition_offsets_[query_vertex][i];
    for (; j < offset; ++j) {
      if (merged_candidates_[query_vertex][j] == INVALID_VERTEX_ID) {
        invalid_sum++;
      } else {
        merged_candidates_[query_vertex][valid_idx++] = merged_candidates_[query_vertex][j];
      }
    }
    CHECK_LE(invalid_sum, offset) << "Error: invalid vertex number greater than offset.";
    CHECK_LE(invalid_sum - last_invalid_sum, per_partition_candidate_cardinality_[i - 1][query_vertex])
        << "Error: per partition invalid vertex number greater than cardinality.";

    offset -= invalid_sum;
    per_partition_candidate_cardinality_[i - 1][query_vertex] -= invalid_sum - last_invalid_sum;
    merged_candidates_[query_vertex].resize(valid_idx);
    last_invalid_sum = invalid_sum;
  }
}

}  // namespace circinus

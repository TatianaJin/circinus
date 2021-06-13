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

std::unique_ptr<Result> Result::newExecutionResult() { return std::make_unique<ExecutionResult>(); }

void CandidateResult::collect(TaskBase* task) {
  auto scan = dynamic_cast<ScanTask*>(task);
  if (scan == nullptr) {
    return;
  }
  auto& shard_candidates = scan->getScanContext().candidates;
  if (!shard_candidates.empty()) {
    candidates_[task->getTaskId()].push_back(std::move(shard_candidates));
  }
}

void CandidateResult::merge() {
  merged_candidates_.resize(candidates_.size());
  for (uint32_t i = 0; i < candidates_.size(); ++i) {
    for (uint32_t j = 0; j < candidates_[i].size(); ++j) {
      merged_candidates_[i].insert(merged_candidates_[i].end(), candidates_[i][j].begin(), candidates_[i][j].end());
    }
  }
}

void CandidateResult::remove_invalid(QueryVertexID query_vertex) {
  merged_candidates_[query_vertex].erase(
      std::remove_if(merged_candidates_[query_vertex].begin(), merged_candidates_[query_vertex].end(),
                     [invalid = INVALID_VERTEX_ID](VertexID vid) { return vid == invalid; }),
      merged_candidates_[query_vertex].end());
}

void ExecutionResult::collect(TaskBase* task) {
  // FIXME(tatiana)
}

}  // namespace circinus

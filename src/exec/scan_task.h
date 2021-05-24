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
#include "graph/graph.h"
#include "ops/scans.h"

namespace circinus {

class ScanTask : public TaskBase {
 private:
  ScanContext scan_context_;
  const Scan* scan_;  // not owned
  const Graph* graph_;

 public:
  ScanTask(QueryId query_id, TaskId task_id, uint32_t shard_id, const Scan* scan, const Graph* graph)
      : TaskBase(query_id, task_id), scan_context_(scan->initScanContext(shard_id)), scan_(scan), graph_(graph) {}

  ScanContext& getScanContext() { return scan_context_; }

  void run() override { scan_->scan(graph_, &scan_context_); }
};

}  // namespace circinus

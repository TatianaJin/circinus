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
#include "graph/graph_partition.h"
#include "ops/scans.h"

namespace circinus {

class ScanTask : public TaskBase {
 private:
  ScanContext scan_context_;
  const Scan* scan_;  // not owned
  const GraphBase* graph_;
  uint32_t partition_ = ~0u;

 public:
  ScanTask(QueryId query_id, TaskId task_id, uint32_t shard_id, const Scan* scan, const GraphBase* graph)
      : TaskBase(query_id, task_id), scan_context_(scan->initScanContext(shard_id)), scan_(scan), graph_(graph) {
    CHECK(dynamic_cast<const Graph*>(graph_) != nullptr);
    DLOG(INFO) << "ScanTask " << task_id << '.' << shard_id << " [" << scan_context_.scan_offset << ", "
               << scan_context_.scan_end << ")";
  }

  ScanTask(QueryId query_id, TaskId task_id, uint32_t shard_id, const Scan* scan, const GraphBase* graph,
           uint32_t partition)
      : TaskBase(query_id, task_id),
        scan_context_(scan->initScanContext(shard_id)),
        scan_(scan),
        graph_(graph),
        partition_(partition) {
    CHECK(dynamic_cast<const ReorderedPartitionedGraph*>(graph_) != nullptr);
  }

  inline ScanContext& getScanContext() { return scan_context_; }
  inline uint32_t getPartition() const { return partition_; }

  const GraphBase* getDataGraph() const override { return graph_; }

  void run(uint32_t executor_idx) override {
    if (partition_ == ~0u) {
      scan_->scan(dynamic_cast<const Graph*>(graph_), &scan_context_);
    } else {
      // zero-copy partition view
      GraphPartition graph_partition(dynamic_cast<const ReorderedPartitionedGraph*>(graph_), partition_);
      scan_->scan(&graph_partition, &scan_context_);
    }
  }
};

}  // namespace circinus

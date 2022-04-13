#pragma once

#include <utility>

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
  ScanTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
           uint32_t shard_id, const Scan* scan, const GraphBase* graph, std::pair<QueryVertexID, VertexID> seed)
      : TaskBase(query_id, task_id, stop_time),
        scan_context_(scan->initScanContext(task_id, shard_id, seed)),
        scan_(scan),
        graph_(graph) {
    CHECK(dynamic_cast<const Graph*>(graph_) != nullptr);
    DLOG(INFO) << "ScanTask " << task_id << '.' << shard_id << " [" << scan_context_.scan_offset << ", "
               << scan_context_.scan_end << ")";
  }

  ScanTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
           uint32_t shard_id, const Scan* scan, const GraphBase* graph, uint32_t partition,
           std::pair<QueryVertexID, VertexID> seed)
      : TaskBase(query_id, task_id, stop_time),
        scan_context_(scan->initScanContext(task_id, shard_id, seed)),
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

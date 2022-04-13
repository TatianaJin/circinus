#pragma once

#include <vector>

#include "exec/task.h"
#include "graph/graph.h"
#include "ops/filters.h"

namespace circinus {

class NeighborhoodFilterTask : public TaskBase {
 private:
  FilterContext filter_context_;
  const NeighborhoodFilter* filter_;  // not owned
  const GraphBase* graph_;
  std::vector<std::vector<VertexID>>* candidates_;

 public:
  NeighborhoodFilterTask(QueryId query_id, TaskId task_id, std::chrono::time_point<std::chrono::steady_clock> stop_time,
                         uint32_t shard_id, const NeighborhoodFilter* filter, const GraphBase* graph,
                         std::vector<std::vector<VertexID>>* candidates)
      : TaskBase(query_id, task_id, stop_time),
        filter_context_(filter->initFilterContext(shard_id)),
        filter_(filter),
        graph_(graph),
        candidates_(candidates) {}

  FilterContext& getFilterContext() { return filter_context_; }

  const NeighborhoodFilter* getFilter() const { return filter_; }

  const GraphBase* getDataGraph() const override { return graph_; }

  void run(uint32_t executor_idx) override { filter_->filter(graph_, candidates_, &filter_context_); }
};

}  // namespace circinus

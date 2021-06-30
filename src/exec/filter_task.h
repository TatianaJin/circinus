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
  NeighborhoodFilterTask(QueryId query_id, TaskId task_id, uint32_t shard_id, const NeighborhoodFilter* filter,
                         const GraphBase* graph, std::vector<std::vector<VertexID>>* candidates)
      : TaskBase(query_id, task_id),
        filter_context_(filter->initFilterContext(shard_id)),
        filter_(filter),
        graph_(graph),
        candidates_(candidates) {}

  FilterContext& getFilterContext() { return filter_context_; }

  const NeighborhoodFilter* getFilter() const { return filter_; }

  const GraphBase* getDataGraph() const override { return graph_; }

  void run() override { filter_->filter(graph_, candidates_, &filter_context_); }
};

}  // namespace circinus

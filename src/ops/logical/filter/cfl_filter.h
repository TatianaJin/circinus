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

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/logical/filter/filter.h"

namespace circinus {

class NeighborhoodFilter;  // forward declaration

class LogicalCFLFilter : public LogicalNeighborhoodFilter {
 private:
  QueryVertexID start_vertex_;
  std::vector<TreeNode> bfs_tree_;
  std::vector<QueryVertexID> bfs_order_;
  std::vector<QueryVertexID> level_offset_;
  uint32_t level_num_;

 public:
  LogicalCFLFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                   const std::vector<VertexID>& candidate_size);

  std::vector<std::unique_ptr<NeighborhoodFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                       ExecutionConfig& exec) override;
};

}  // namespace circinus

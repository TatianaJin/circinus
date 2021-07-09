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

#include <memory>
#include <unordered_map>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "ops/logical/filter/filter.h"

namespace circinus {

class NeighborhoodFilter;  // forward declaration

class LogicalTSOFilter : public LogicalNeighborhoodFilter {
 private:
  QueryVertexID start_vertex_;
  std::vector<TreeNode> tree_;
  std::vector<QueryVertexID> dfs_order_;

 public:
  LogicalTSOFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                   const std::vector<VertexID>& candidate_size);

  QueryVertexID getStartVertex(const GraphMetadata& metadata, const QueryGraph* query_graph_,
                               const std::vector<VertexID>& candidate_size);

  const std::vector<QueryVertexID>& getDfsOrder() const { return dfs_order_; }

  const std::vector<TreeNode>& getTree() const { return tree_; }

  std::vector<std::unique_ptr<NeighborhoodFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                       ExecutionConfig& exec) override;
};

}  // namespace circinus

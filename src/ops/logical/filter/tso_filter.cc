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

#include "ops/logical/filter/tso_filter.h"

#include <memory>

#include "ops/filters/filter.h"
#include "ops/order/tso_order.h"
#include "utils/utils.h"

namespace circinus {

LogicalTSOFilter::LogicalTSOFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                                   const std::vector<VertexID>& candidate_size)
    : LogicalNeighborhoodFilter(query_graph) {
  TSOOrder tso_order;
  start_vertex_ = tso_order.getStartVertex(metadata, query_graph, candidate_size);
  std::vector<QueryVertexID> bfs_order_;
  bfs(query_graph, start_vertex_, tree_, bfs_order_);
  dfs(start_vertex_, tree_, dfs_order_);
}

std::vector<std::unique_ptr<NeighborhoodFilter>> LogicalTSOFilter::toPhysicalOperators(const GraphMetadata& metadata,
                                                                                       ExecutionConfig& exec) {
  std::vector<std::unique_ptr<NeighborhoodFilter>> ret;

  for (QueryVertexID query_vertex : dfs_order_) {
    if (start_vertex_ == query_vertex) {
      continue;
    }
    TreeNode& node = tree_[query_vertex];
    LOG(INFO) << query_vertex << " " << node.parent_;
    ret.emplace_back(std::make_unique<NeighborhoodFilter>(exec, query_graph_, query_vertex, node.parent_));
  }

  for (auto it = dfs_order_.rbegin(); it != dfs_order_.rend(); ++it) {
    TreeNode& node = tree_[*it];
    if (node.children_.size() > 0) {
      ret.emplace_back(std::make_unique<NeighborhoodFilter>(exec, query_graph_, *it, node.children_));
    }
  }

  return ret;
}

}  // namespace circinus

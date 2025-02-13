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

#include "ops/logical/filter/dpiso_filter.h"

#include <memory>

#include "exec/execution_config.h"
#include "ops/filters/filter.h"
#include "ops/order/dpiso_order.h"
#include "utils/utils.h"

namespace circinus {

LogicalDPISOFilter::LogicalDPISOFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                                       const std::vector<VertexID>& candidate_size)
    : LogicalNeighborhoodFilter(query_graph) {
  OrderBase dpiso_order;
  start_vertex_ = dpiso_order.getStartVertex(metadata, query_graph, candidate_size);
  uint32_t query_vertices_num = query_graph->getNumVertices();
  bfs(query_graph, start_vertex_, tree_, bfs_order_);
  std::vector<uint32_t> order_index(query_vertices_num);
  for (uint32_t i = 0; i < query_vertices_num; ++i) {
    QueryVertexID query_vertex = bfs_order_[i];
    order_index[query_vertex] = i;
  }

  for (uint32_t i = 0; i < query_vertices_num; ++i) {
    QueryVertexID u = bfs_order_[i];

    const auto& u_nbrs = query_graph->getOutNeighbors(u);
    for (uint32_t j = 0; j < u_nbrs.second; ++j) {
      QueryVertexID u_nbr = u_nbrs.first[j];
      if (order_index[u_nbr] < order_index[u]) {
        tree_[u].bn_.emplace_back(u_nbr);
      } else {
        tree_[u].fn_.emplace_back(u_nbr);
      }
    }
  }
}

std::vector<std::unique_ptr<NeighborhoodFilter>> LogicalDPISOFilter::toPhysicalOperators(const GraphMetadata& metadata,
                                                                                         ExecutionConfig& exec) {
  std::vector<std::unique_ptr<NeighborhoodFilter>> ret;
  for (uint32_t refine_time = 0; refine_time < 3; ++refine_time) {
    if (refine_time % 2 == 0) {
      for (QueryVertexID query_vertex : bfs_order_) {
        TreeNode& node = tree_[query_vertex];
        if (node.bn_.size() > 0) {
          ret.emplace_back(std::make_unique<NeighborhoodFilter>(exec, query_graph_, query_vertex, node.bn_));
          // merge output operator
        }
      }
    } else {
      for (auto it = bfs_order_.rbegin(); it != bfs_order_.rend(); ++it) {
        TreeNode& node = tree_[*it];
        if (node.fn_.size() > 0) {
          ret.emplace_back(std::make_unique<NeighborhoodFilter>(exec, query_graph_, *it, node.fn_));
          // merge output operator
        }
      }
    }
  }
  return ret;
}

}  // namespace circinus

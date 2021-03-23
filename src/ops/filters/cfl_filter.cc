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

#include <algorithm>
#include <vector>

#include "ops/filters/cfl_filter.h"
#include "utils/utils.h"

namespace circinus {

CFLFilter::CFLFilter(const QueryGraph* query_graph, const Graph* data_graph, QueryVertexID start_vertex)
    : FilterBase(query_graph, data_graph) {
  start_vertex_ = start_vertex;
  uint32_t query_vertices_num = query_graph->getNumVertices();
  bfs_tree_.resize(query_vertices_num);
  bfs_order_.reserve(query_vertices_num);
  bfs(query_graph, start_vertex, bfs_tree_, bfs_order_);

  std::vector<uint32_t> order_index(query_vertices_num);
  for (uint32_t i = 0; i < query_vertices_num; ++i) {
    QueryVertexID query_vertex = bfs_order_[i];
    order_index[query_vertex] = i;
  }

  level_num_ = -1;
  level_offset_.reserve(query_vertices_num + 1);

  for (uint32_t i = 0; i < query_vertices_num; ++i) {
    QueryVertexID u = bfs_order_[i];

    if (bfs_tree_[u].level_ != level_num_) {
      level_num_ += 1;
      level_offset_[level_num_] = 0;
    }

    level_offset_[level_num_] += 1;

    const auto& u_nbrs = query_graph->getOutNeighbors(u);
    for (uint32_t j = 0; j < u_nbrs.second; ++j) {
      QueryVertexID u_nbr = u_nbrs.first[j];

      if (bfs_tree_[u].level_ == bfs_tree_[u_nbr].level_) {
        if (order_index[u_nbr] < order_index[u]) {
          bfs_tree_[u].bn_.emplace_back(u_nbr);
        } else {
          bfs_tree_[u].fn_.emplace_back(u_nbr);
        }
      } else if (bfs_tree_[u].level_ > bfs_tree_[u_nbr].level_) {
        bfs_tree_[u].bn_.emplace_back(u_nbr);
      } else {
        bfs_tree_[u].under_level_.emplace_back(u_nbr);
      }
    }
  }

  level_num_ += 1;

  for (uint32_t i = level_num_; i >= 1; --i) {
    level_offset_[i] = level_offset_[i - 1];
  }
  level_offset_[0] = 0;
  for (uint32_t i = 1; i <= level_num_; ++i) {
    level_offset_[i] += level_offset_[i - 1];
  }
}

void CFLFilter::Filter(std::vector<std::vector<VertexID>>& candidates) {
  for (uint32_t i = 0; i < level_num_; ++i) {
    for (uint32_t j = level_offset_[i]; j < level_offset_[i + 1]; ++j) {
      QueryVertexID query_vertex = bfs_order_[j];
      TreeNode& node = bfs_tree_[query_vertex];
      if (node.bn_.size() > 0) {
        pruneByPivotVertices(query_vertex, node.bn_, candidates);
      }
    }

    for (int j = (int)level_offset_[i + 1] - 1; j >= (int)level_offset_[i]; --j) {
      QueryVertexID query_vertex = bfs_order_[j];
      TreeNode& node = bfs_tree_[query_vertex];
      if (node.fn_.size() > 0) {
        pruneByPivotVertices(query_vertex, node.fn_, candidates);
      }
    }
  }

  for (int i = (int)level_num_ - 2; i >= 0; --i) {
    for (uint32_t j = level_offset_[i]; j < level_offset_[i + 1]; ++j) {
      QueryVertexID query_vertex = bfs_order_[j];
      TreeNode& node = bfs_tree_[query_vertex];
      if (node.under_level_.size() > 0) {
        pruneByPivotVertices(query_vertex, node.under_level_, candidates);
      }
    }
  }
}

}  // namespace circinus

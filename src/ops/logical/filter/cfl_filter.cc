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
#include <memory>
#include <vector>

#include "exec/execution_config.h"
#include "ops/filters/filter.h"
#include "ops/logical/filter/cfl_filter.h"
#include "ops/order/cfl_order.h"
#include "utils/utils.h"

namespace circinus {

LogicalCFLFilter::LogicalCFLFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                                   const std::vector<VertexID>& candidate_size)
    : LogicalNeighborhoodFilter(query_graph),two_core_solver_(query_graph) {
  CFLOrder cfl_order;
  start_vertex_ = cfl_order.getStartVertex(metadata, query_graph, candidate_size);
  LOG(INFO) << "start_vertex " << start_vertex_;
  uint32_t query_vertices_num = query_graph->getNumVertices();
  bfs_tree_.resize(query_vertices_num);
  bfs_order_.reserve(query_vertices_num);
  bfs(query_graph, start_vertex_, bfs_tree_, bfs_order_);

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

std::vector<std::unique_ptr<NeighborhoodFilter>> LogicalCFLFilter::toPhysicalOperators(const GraphMetadata& metadata,
                                                                                       ExecutionConfig& exec) {
  std::vector<std::unique_ptr<NeighborhoodFilter>> ret;
  for (uint32_t i = 0; i < level_num_; ++i) {
    for (uint32_t j = level_offset_[i]; j < level_offset_[i + 1]; ++j) {
      QueryVertexID query_vertex = bfs_order_[j];
      TreeNode& node = bfs_tree_[query_vertex];
      if (node.bn_.size() > 0) {
        ret.emplace_back(std::make_unique<NeighborhoodFilter>(exec, query_graph_, query_vertex, node.bn_));
        // merge output operator
      }
    }

    for (int j = (int)level_offset_[i + 1] - 1; j >= (int)level_offset_[i]; --j) {
      QueryVertexID query_vertex = bfs_order_[j];
      TreeNode& node = bfs_tree_[query_vertex];
      if (node.fn_.size() > 0) {
        ret.emplace_back(std::make_unique<NeighborhoodFilter>(exec, query_graph_, query_vertex, node.fn_));
        // merge output operator
      }
    }
  }

  for (int i = (int)level_num_ - 2; i >= 0; --i) {
    for (uint32_t j = level_offset_[i]; j < level_offset_[i + 1]; ++j) {
      QueryVertexID query_vertex = bfs_order_[j];
      TreeNode& node = bfs_tree_[query_vertex];
      if (node.under_level_.size() > 0) {
        ret.emplace_back(std::make_unique<NeighborhoodFilter>(exec, query_graph_, query_vertex, node.under_level_));
        // merge output operator
      }
    }
  }

  return ret;
}
std::vector<QueryVertexID> LogicalCFLFilter::getTopThree(const GraphMetadata& metadata, const QueryGraph* q) {
  auto rank_compare = [](std::pair<QueryVertexID, double> l, std::pair<QueryVertexID, double> r) {
    return l.second < r.second;
  };
  std::priority_queue<std::pair<QueryVertexID, double>, std::vector<std::pair<QueryVertexID, double>>,
                      decltype(rank_compare)>
      rank_queue(rank_compare);
  const auto& core_table = two_core_solver_.get2CoreTable();
  uint32_t core_size = two_core_solver_.getCoreSize();

  for (QueryVertexID query_vertex = 0; query_vertex < q->getNumVertices(); ++query_vertex) {
    if (core_size == 0 || two_core_solver_.isInCore(core_table, query_vertex)) {
      double degree = q->getVertexOutDegree(query_vertex);
      uint32_t frequency = metadata.getLabelFrequency(q->getVertexLabel(query_vertex));
      double rank = frequency / degree;
      rank_queue.emplace(query_vertex, rank);
    }
  }
  std::vector<QueryVertexID> ret(3);
  // Keep the top-3.
  while (rank_queue.size() > 3) {
    rank_queue.pop();
  }
  for (uint32_t i = 0; i < 3; ++i) {
    ret[i] = rank_queue.top().first;
    rank_queue.pop();
  }
  return ret;
}

QueryVertexID LogicalCFLFilter::getStartVertex(const std::vector<QueryVertexID>& query_vertices,
                                       const std::vector<VertexID>& candidate_size, const QueryGraph& query_graph,
                                       const GraphMetadata& metadata) {
  QueryVertexID start_vertex = 0;
  double min_score = metadata.getNumVertices() + 1;
  for (auto query_vertex : query_vertices) {
    double cur_score = candidate_size[query_vertex] / (double)query_graph.getVertexOutDegree(query_vertex);
    if (cur_score < min_score) {
      start_vertex = query_vertex;
      min_score = cur_score;
    }
  }
  return start_vertex;
}

QueryVertexID LogicalCFLFilter::getStartVertex(const GraphMetadata& metadata, const QueryGraph* query_graph,
                                       const std::vector<VertexID>& candidate_size) {
  auto root_candidates = getTopThree(metadata, query_graph);
  return getStartVertex(root_candidates, {candidate_size.begin(), candidate_size.end()}, *query_graph, metadata);
}

}  // namespace circinus

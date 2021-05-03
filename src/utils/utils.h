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

#include <cstring>
#include <queue>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "graph/types.h"

namespace circinus {

const uint64_t INVALID_VERTEX_ID = 0xffffffffffffffff;

using circinus::QueryGraph;

static void dfs(QueryVertexID cur_vertex, const std::vector<TreeNode>& dfs_tree,
                std::vector<QueryVertexID>& dfs_order) {
  dfs_order.emplace_back(cur_vertex);
  for (QueryVertexID v : dfs_tree[cur_vertex].children_) {
    dfs(v, dfs_tree, dfs_order);
  }
}

inline static void bfs(const QueryGraph* query_graph, QueryVertexID start_vertex, std::vector<TreeNode>& bfs_tree,
                       std::vector<QueryVertexID>& bfs_order) {
  QueryVertexID vertex_num = query_graph->getNumVertices();
  bfs_tree.resize(vertex_num);
  bfs_order.resize(vertex_num);
  std::queue<QueryVertexID> bfs_queue;
  std::vector<bool> visited(vertex_num, false);

  uint32_t visited_vertex_count = 0;
  bfs_queue.push(start_vertex);
  visited[start_vertex] = true;
  bfs_tree[start_vertex].level_ = 0;

  while (!bfs_queue.empty()) {
    const QueryVertexID u = bfs_queue.front();
    bfs_queue.pop();
    bfs_order[visited_vertex_count++] = u;

    const auto& u_nbrs = query_graph->getOutNeighbors(u);
    for (uint32_t i = 0; i < u_nbrs.second; ++i) {
      QueryVertexID u_nbr = u_nbrs.first[i];

      if (!visited[u_nbr]) {
        bfs_queue.push(u_nbr);
        visited[u_nbr] = true;
        bfs_tree[u_nbr].parent_ = u;
        bfs_tree[u_nbr].level_ = bfs_tree[u].level_ + 1;
        bfs_tree[u].children_.emplace_back(u_nbr);
      }
    }
  }
}

inline static bool semiperfectBipartiteMatching(uint32_t* col_ptrs, const std::vector<uint32_t>& col_ids, uint32_t n,
                                                uint32_t m) {
  std::vector<int> match(n, -1);
  std::vector<int> row_match(m, -1);
  std::vector<int> visited(m, -1);
  std::vector<int> previous(m, 0);
  std::queue<uint32_t> q;

  for (uint32_t i = 0; i < n; ++i) {
    if (match[i] == -1 && col_ptrs[i] != col_ptrs[i + 1]) {
      // clear queue
      while (!q.empty()) q.pop();
      q.push(i);
      previous[i] = -1;

      bool flag = false;
      while (!q.empty() && !flag) {
        uint32_t u = q.front();
        q.pop();
        for (uint32_t ptr = col_ptrs[u]; ptr < col_ptrs[u + 1] && !flag; ptr++) {
          uint32_t v = col_ids[ptr];
          if (visited[v] == -1 || ((uint32_t)visited[v]) != i) {
            visited[v] = i;
            if (row_match[v] >= 0) {
              q.push(row_match[v]);
              previous[row_match[v]] = u;
            } else {
              flag = true;
              int d = u, e = v;
              while (d != -1) {
                int t = match[d];
                match[d] = e;
                row_match[e] = d;
                d = previous[d];
                e = t;
              }
            }
          }
        }
      }
    }
    if (match[i] == -1) {
      return false;
    }
  }
  return true;
}

static inline uint64_t getNumSubgraphs(const std::vector<CompressedSubgraphs>& groups, uint32_t start, uint32_t end) {
  uint64_t count = 0;
  for (uint32_t i = start; i < end; ++i) {
    count += groups[i].getNumSubgraphs();
  }
  return count;
}

}  // namespace circinus

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

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "graph/types.h"

namespace circinus {

const uint32_t INVALID_VERTEX_ID = 0xffffffff;

using circinus::QueryGraph;

static void dfs(QueryVertexID cur_vertex, const std::vector<TreeNode>& dfs_tree,
                std::vector<QueryVertexID>& dfs_order) {
  dfs_order.emplace_back(cur_vertex);
  for (QueryVertexID v : dfs_tree[cur_vertex].children_) {
    dfs(v, dfs_tree, dfs_order);
  }
}

static void bfs(const QueryGraph* query_graph, QueryVertexID start_vertex, std::vector<TreeNode>& bfs_tree,
                std::vector<QueryVertexID>& bfs_order) {
  QueryVertexID vertex_num = query_graph->getNumVertices();

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

static void precompute(uint32_t* col_ptrs, const std::vector<uint32_t>& col_ids, std::vector<int>& match,
                       std::vector<int>& row_match, uint32_t n, uint32_t m) {
  for (uint32_t i = 0; i < n; i++) {
    uint32_t s_ptr = col_ptrs[i];
    uint32_t e_ptr = col_ptrs[i + 1];
    for (uint32_t ptr = s_ptr; ptr < e_ptr; ptr++) {
      uint32_t r_id = col_ids[ptr];
      if (row_match[r_id] == -1) {
        match[i] = r_id;
        row_match[r_id] = i;
        break;
      }
    }
  }
}

static bool semiperfectBipartiteMatching(uint32_t* col_ptrs, const std::vector<uint32_t>& col_ids, uint32_t n,
                                         uint32_t m) {
  std::vector<int> match(n, -1);
  std::vector<int> row_match(m, -1);
  std::vector<int> visited(m, -1);
  std::vector<int> previous(m, 0);
  std::queue<uint32_t> q;
  precompute(col_ptrs, col_ids, match, row_match, n, m);

  for (uint32_t i = 0; i < n; i++) {
    if (match[i] == -1 && col_ptrs[i] != col_ptrs[i + 1]) {
      // clear queue
      while (!q.empty()) q.pop();
      q.push(i);

      while (!q.empty()) {
        uint32_t queue_col = q.front();
        q.pop();
        uint32_t eptr = col_ptrs[queue_col + 1];
        for (uint32_t ptr = col_ptrs[queue_col]; ptr < eptr; ptr++) {
          int row = col_ids[ptr];

          if (visited[row] != (int)i) {
            previous[row] = queue_col;
            visited[row] = i;

            int col = row_match[row];

            if (col == -1) {
              // Find an augmenting path. Then, trace back and modify the augmenting path.
              while (row != -1) {
                col = previous[row];
                int temp = match[col];
                match[col] = row;
                row_match[row] = col;
                row = temp;
              }
              break;
            } else {
              // Continue to construct the match.
              q.push(col);
            }
          }
        }
      }

      if (match[i] == -1) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace circinus

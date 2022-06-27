#pragma once

#include <cstring>
#include <memory>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "graph/types.h"

namespace circinus {

#define toSeconds(start, end) \
  (((double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1e9)

#define toMilliseconds(start, end) \
  (((double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1e6)

const VertexID INVALID_VERTEX_ID = std::numeric_limits<VertexID>::max();

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

template <typename V>
inline std::string toString(const std::vector<V>& vec) {
  std::stringstream ss;
  for (auto v : vec) {
    ss << ' ' << v;
  }
  return ss.str();
}

template <typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& ss, const std::pair<T1, T2>& p) {
  return ss << p.first << '-' << p.second;
}

template <typename ForwardIt>
inline std::string toString(ForwardIt begin, ForwardIt end) {
  std::stringstream ss;
  for (auto iter = begin; iter != end; ++iter) {
    ss << ' ' << *iter;
  }
  return ss.str();
}

#ifdef __GNUG__
#include <cxxabi.h>
#include <cstdlib>
template <typename T>
std::string getTypename(const T& t) {
  auto name = typeid(t).name();
  int status = 0;
  // enable c++11 by passing the flag -std=c++11 to g++
  std::unique_ptr<char, void (*)(void*)> res{abi::__cxa_demangle(name, NULL, NULL, &status), std::free};
  return (status == 0) ? res.get() : name;
}
#else
template <typename T>
std::string getTypename(const T& t) {
  return typeid(t).name();
}
#endif

}  // namespace circinus

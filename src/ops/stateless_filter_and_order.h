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
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#ifdef WITH_GPERF
#include "gperftools/profiler.h"
#endif
#include "gtest/gtest.h"

#include "exec/thread_pool.h"
#include "graph/bipartite_graph.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters.h"
#include "ops/logical_filters.h"
#include "ops/operators.h"
#include "ops/scans.h"
#include "ops/types.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/hashmap.h"
#include "utils/profiler.h"
#include "utils/utils.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionConfig;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::BipartiteGraph;
using circinus::GraphType;
using circinus::GraphMetadata;
using circinus::NaivePlanner;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::Task;
using circinus::ThreadPool;
using circinus::VertexID;
using circinus::Profiler;
using circinus::CoverNode;
using circinus::QueryType;
using circinus::TraverseOperator;
using circinus::INVALID_VERTEX_ID;

// logical filter
using circinus::LogicalCFLFilter;
using circinus::LogicalGQLFilter;
using circinus::LogicalNLFFilter;
using circinus::LogicalTSOFilter;
using circinus::LogicalDPISOFilter;
using circinus::LogicalNeighborhoodFilter;

// physical filter
using circinus::NeighborhoodFilter;
using circinus::NLFFilter;
using circinus::GQLFilter;
namespace circinus {
class StatelessFilterAndOrder {
 public:
  static QueryVertexID selectGQLStartVertex(const QueryGraph* query_graph,
                                            const std::vector<VertexID>& candidates_cardinality) {
    QueryVertexID start_vertex = 0;
    QueryVertexID qg_v_cnt = query_graph->getNumVertices();
    for (QueryVertexID i = 1; i < qg_v_cnt; ++i) {
      QueryVertexID cur_vertex = i;
      VertexID size1 = candidates_cardinality[cur_vertex], size2 = candidates_cardinality[start_vertex];
      if (size1 < size2) {
        start_vertex = cur_vertex;
      } else if (size1 == size2 &&
                 query_graph->getVertexOutDegree(cur_vertex) > query_graph->getVertexOutDegree(start_vertex)) {
        start_vertex = cur_vertex;
      }
    }
    return start_vertex;
  }

  static void updateValidVertices(const QueryGraph* query_graph, QueryVertexID query_vertex, std::vector<bool>& visited,
                                  std::vector<bool>& adjacent) {
    visited[query_vertex] = true;
    auto[nbrs, cnt] = query_graph->getOutNeighbors(query_vertex);
    for (uint32_t i = 0; i < cnt; ++i) {
      adjacent[nbrs[i]] = true;
    }
  }

  static std::vector<QueryVertexID> getGQLOrder(const Graph* data_graph, const QueryGraph* query_graph,
                                                const std::vector<VertexID>& candidates_cardinality) {
    QueryVertexID qg_v_cnt = query_graph->getNumVertices();
    std::vector<bool> visited_vertices(qg_v_cnt, false);
    std::vector<bool> adjacent_vertices(qg_v_cnt, false);
    std::vector<QueryVertexID> order(qg_v_cnt);

    QueryVertexID start_vertex = selectGQLStartVertex(query_graph, candidates_cardinality);
    order[0] = start_vertex;
    updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

    for (QueryVertexID i = 1; i < qg_v_cnt; ++i) {
      QueryVertexID next_vertex;
      QueryVertexID min_value = data_graph->getNumVertices() + 1;
      for (QueryVertexID j = 0; j < qg_v_cnt; ++j) {
        QueryVertexID cur_vertex = j;
        if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
          size_t cnt = candidates_cardinality[cur_vertex];
          if (cnt < min_value) {
            min_value = cnt;
            next_vertex = cur_vertex;
          } else if (cnt == min_value &&
                     query_graph->getVertexOutDegree(cur_vertex) > query_graph->getVertexOutDegree(next_vertex)) {
            next_vertex = cur_vertex;
          }
        }
      }
      updateValidVertices(next_vertex, visited_vertices, adjacent_vertices);
      order[i] = next_vertex;
    }
    return order;
  }

  static std::vector<QueryVertexID> getDPISOOrder(const QueryGraph* query_graph,
                                                  const LogicalDPISOFilter& logical_filter) {
    const std::vector<QueryVertexID>& bfs_order = logical_filter.getBfsOrder();
    QueryVertexID qg_v_cnt = query_graph->getNumVertices();
    std::vector<QueryVertexID> order(qg_v_cnt);
    for (QueryVertexID i = 0; i < qg_v_cnt; ++i) {
      order[i] = bfs_order[i];
    }
    return order;
  }

  static void estimatePathEmbeddsingsNum(std::vector<QueryVertexID>& path,
                                         std::vector<size_t>& estimated_embeddings_num,
                                         std::vector<std::vector<VertexID>> candidates_sets,
                                         std::map<std::pair<QueryVertexID, QueryVertexID>, BipartiteGraph> bg_map) {
    assert(path.size() > 1);
    std::vector<size_t> parent;
    std::vector<size_t> children;

    size_t begin = path.size() - 2, end = path.size() - 1;

    estimated_embeddings_num.resize(path.size() - 1);
    auto last_edge = bg_map[{path[begin], path[end]}];
    children.resize(last_edge->getNumVertices());

    size_t sum = 0;
    for (auto& v : candidates_sets[path[begin]]) {
      int offset = last_edge->getOffset(v);
      children[offset] = last_edge->getVertexOutDegree(v);
      sum += children[offset];
    }

    estimated_embeddings_num[begin] = sum;

    for (int i = begin; i >= 1; --i) {
      begin = i - 1;
      end = i;
      auto edge = bg_map[{path[begin], path[end]}];
      parent.resize(edge->getNumVertices());

      sum = 0;
      for (auto& v : candidates_sets[path[begin]]) {
        size_t local_sum = 0;
        auto[nbrs, cnt] = edge->getOutNeighbors(v);
        for (uint32_t j = 0; j < cnt; ++j) {
          auto nbr = nbrs[j];
          local_sum += children[last_edge->getOffset(nbr)];
        }
        parent[edge->getOffset(v)] = local_sum;
        sum += local_sum;
      }

      estimated_embeddings_num[i - 1] = sum;
      parent.swap(children);

      last_edge = edge;
    }
  }

  static QueryVertexID generateNoneTreeEdgesCount(const QueryGraph* query_graph, const std::vector<TreeNode>& tree_node,
                                                  std::vector<QueryVertexID>& path) {
    auto non_tree_edge_count = query_graph->getVertexOutDegree(path[0]) - tree_node[path[0]].children_.size();
    for (size_t i = 1; i < path.size(); ++i) {
      auto vertex = path[i];
      non_tree_edge_count += query_graph->getVertexOutDegree(vertex) - tree_node[vertex].children_.size() - 1;
    }

    return non_tree_edge_count;
  }

  static void generateRootToLeafPaths(const std::vector<TreeNode>& tree_node, QueryVertexID cur_vertex,
                                      std::vector<QueryVertexID>& cur_path,
                                      std::vector<std::vector<QueryVertexID>>& paths) {
    auto& cur_node = tree_node[cur_vertex];
    cur_path.push_back(cur_vertex);
    if (cur_node.children_.size() == 0) {
      paths.emplace_back(cur_path);
    } else {
      for (auto child : cur_node.children_) {
        generateRootToLeafPaths(tree_node, child, cur_path, paths);
      }
    }
    cur_path.pop_back();
  }

  static std::vector<QueryVertexID> getTSOOrder(
      const QueryGraph* query_graph, const LogicalTSOFilter& logical_filter,
      std::vector<std::vector<VertexID>> candidates_sets,
      std::map<std::pair<QueryVertexID, QueryVertexID>, BipartiteGraph> bg_map) {
    const std::vector<TreeNode>& tree = logical_filter.getTree();
    const std::vector<QueryVertexID>& dfs_order = logical_filter.getDfsOrder();

    QueryVertexID qg_v_cnt = query_graph->getNumVertices();
    std::vector<std::vector<QueryVertexID>> paths;
    paths.reserve(qg_v_cnt);

    std::vector<QueryVertexID> single_path;
    single_path.reserve(qg_v_cnt);

    generateRootToLeafPaths(tree, dfs_order[0], single_path, paths);
    std::vector<std::pair<double, std::vector<QueryVertexID>*>> path_orders;
    for (auto& path : paths) {
      std::vector<size_t> estimated_embeddings_num;
      QueryVertexID non_tree_edges_count = generateNoneTreeEdgesCount(tree, path);
      estimatePathEmbeddsingsNum(path, estimated_embeddings_num, candidates_sets, bg_map);
      double score = estimated_embeddings_num[0] / (double)(non_tree_edges_count + 1);
      path_orders.emplace_back(std::make_pair(score, &path));
    }
    std::sort(path_orders.begin(), path_orders.end(),
              [](std::pair<double, std::vector<QueryVertexID>*> l, std::pair<double, std::vector<QueryVertexID>*> r) {
                return l.first < r.first;
              });
    std::vector<bool> visited_vertices(qg_v_cnt, false);
    std::vector<QueryVertexID> order;
    order.reserve(qg_v_cnt);
    for (auto& path : path_orders) {
      for (auto v : *(path.second)) {
        if (!visited_vertices[v]) {
          order.push_back(v);
          visited_vertices[v] = true;
        }
      }
    }
    return order;
  }

  static void generateLeaves(const QueryGraph* query_graph, std::vector<QueryVertexID>& leaves) {
    for (QueryVertexID i = 0; i < query_graph->getNumVertices(); ++i) {
      QueryVertexID cur_vertex = i;
      if (query_graph->getVertexOutDegree(cur_vertex) == 1) {
        leaves.push_back(cur_vertex);
      }
    }
  }

  static void generateCorePaths(const QueryGraph* query_graph, const std::vector<TreeNode>& tree_node,
                                QueryVertexID cur_vertex, std::vector<QueryVertexID>& cur_core_path,
                                std::vector<std::vector<QueryVertexID>>& core_paths, const TwoCoreSolver& tcs) {
    const TreeNode& node = tree_node[cur_vertex];
    cur_core_path.push_back(cur_vertex);

    bool is_core_leaf = true;
    for (const auto& child : node.children_) {
      if (tcs.isInCore(child)) {
        generateCorePaths(query_graph, tree_node, child, cur_core_path, core_paths, tcs);
        is_core_leaf = false;
      }
    }

    if (is_core_leaf) {
      core_paths.emplace_back(cur_core_path);
    }
    cur_core_path.pop_back();
  }

  static void generateTreePaths(const QueryGraph* query_graph, const std::vector<TreeNode>& tree_node,
                                QueryVertexID cur_vertex, std::vector<QueryVertexID>& cur_tree_path,
                                std::vector<std::vector<QueryVertexID>>& tree_paths) {
    const TreeNode& node = tree_node[cur_vertex];
    cur_tree_path.push_back(cur_vertex);

    bool is_tree_leaf = true;
    for (auto child : node.children_) {
      if (query_graph->getVertexOutDegree(child) > 1) {
        generateTreePaths(query_graph, tree_node, child, cur_tree_path, tree_paths);
        is_tree_leaf = false;
      }
    }

    if (is_tree_leaf && cur_tree_path.size() > 1) {
      tree_paths.emplace_back(cur_tree_path);
    }
    cur_tree_path.pop_back();
  }

  static std::vector<QueryVertexID> getCFLOrder(
      const Graph* data_graph, const QueryGraph* query_graph, const LogicalCFLFilter& logical_filter,
      std::vector<std::vector<VertexID>> candidates_sets,
      std::map<std::pair<QueryVertexID, QueryVertexID>, BipartiteGraph> bg_map) {
    const std::vector<TreeNode>& tree = logical_filter.getTree();
    const std::vector<QueryVertexID>& bfs_order = logical_filter.getBfsOrder();

    QueryVertexID qg_v_cnt = query_graph->getNumVertices();
    QueryVertexID root_vertex = bfs_order[0];
    std::vector<QueryVertexID> order(qg_v_cnt);
    std::vector<bool> visited_vertices(qg_v_cnt, false);

    std::vector<std::vector<QueryVertexID>> core_paths;
    std::vector<std::vector<std::vector<QueryVertexID>>> forests;
    std::vector<QueryVertexID> leaves;

    generateLeaves(query_graph, leaves);

    const auto& tcs = logical_filter.getTwoCoreSolver();
    const auto& core_table = tcs.get2CoreTable();
    if (tcs.isInCore(root_vertex)) {
      std::vector<QueryVertexID> temp_core_path;
      generateCorePaths(query_graph, tree, root_vertex, temp_core_path, core_paths, tcs);
      for (QueryVertexID i = 0; i < qg_v_cnt; ++i) {
        QueryVertexID cur_vertex = i;
        if (tcs.isInCore(cur_vertex)) {
          std::vector<std::vector<QueryVertexID>> temp_tree_paths;
          std::vector<QueryVertexID> temp_tree_path;
          generateTreePaths(query_graph, tree, cur_vertex, temp_tree_path, temp_tree_paths);
          if (!temp_tree_paths.empty()) {
            forests.emplace_back(temp_tree_paths);
          }
        }
      }
    } else {
      std::vector<std::vector<QueryVertexID>> temp_tree_paths;
      std::vector<QueryVertexID> temp_tree_path;
      generateTreePaths(query_graph, tree, root_vertex, temp_tree_path, temp_tree_paths);
      if (!temp_tree_paths.empty()) {
        forests.emplace_back(temp_tree_paths);
      }
    }

    // Order core paths.
    QueryVertexID selected_vertices_count = 0;
    order[selected_vertices_count++] = root_vertex;
    visited_vertices[root_vertex] = true;

    if (!core_paths.empty()) {
      std::vector<std::vector<size_t>> paths_embededdings_num;
      std::vector<QueryVertexID> paths_non_tree_edge_num;
      for (auto& path : core_paths) {
        QueryVertexID non_tree_edge_num = generateNoneTreeEdgesCount(tree, path);
        paths_non_tree_edge_num.push_back(non_tree_edge_num + 1);

        std::vector<size_t> path_embeddings_num;
        estimatePathEmbeddsingsNum(path, path_embeddings_num, candidates_sets, bg_map);
        paths_embededdings_num.emplace_back(path_embeddings_num);
      }

      // Select the start path.
      double min_value = std::numeric_limits<double>::max();
      size_t selected_path_index = 0;

      for (size_t i = 0; i < core_paths.size(); ++i) {
        double cur_value = paths_embededdings_num[i][0] / (double)paths_non_tree_edge_num[i];

        if (cur_value < min_value) {
          min_value = cur_value;
          selected_path_index = i;
        }
      }

      for (QueryVertexID i = 1; i < core_paths[selected_path_index].size(); ++i) {
        order[selected_vertices_count] = core_paths[selected_path_index][i];
        selected_vertices_count += 1;
        visited_vertices[core_paths[selected_path_index][i]] = true;
      }

      core_paths.erase(core_paths.begin() + selected_path_index);
      paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
      paths_non_tree_edge_num.erase(paths_non_tree_edge_num.begin() + selected_path_index);

      while (!core_paths.empty()) {
        min_value = std::numeric_limits<double>::max();
        selected_path_index = 0;

        for (size_t i = 0; i < core_paths.size(); ++i) {
          QueryVertexID path_root_vertex_idx = 0;
          for (QueryVertexID j = 0; j < core_paths[i].size(); ++j) {
            QueryVertexID cur_vertex = core_paths[i][j];

            if (visited_vertices[cur_vertex]) continue;

            path_root_vertex_idx = j - 1;
            break;
          }

          double cur_value = paths_embededdings_num[i][path_root_vertex_idx] /
                             (double)(candidates_sets_[core_paths[i][path_root_vertex_idx]].size());
          if (cur_value < min_value) {
            min_value = cur_value;
            selected_path_index = i;
          }
        }

        for (QueryVertexID i = 1; i < core_paths[selected_path_index].size(); ++i) {
          if (visited_vertices[core_paths[selected_path_index][i]]) continue;

          order[selected_vertices_count] = core_paths[selected_path_index][i];
          selected_vertices_count += 1;
          visited_vertices[core_paths[selected_path_index][i]] = true;
        }

        core_paths.erase(core_paths.begin() + selected_path_index);
        paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
      }
    }

    // Order tree paths.
    for (auto& tree_paths : forests) {
      std::vector<std::vector<size_t>> paths_embededdings_num;
      for (auto& path : tree_paths) {
        std::vector<size_t> path_embeddings_num;
        estimatePathEmbeddsingsNum(path, path_embeddings_num, candidates_sets, bg_map);
        paths_embededdings_num.emplace_back(path_embeddings_num);
      }

      while (!tree_paths.empty()) {
        double min_value = std::numeric_limits<double>::max();
        QueryVertexID selected_path_index = 0;

        for (size_t i = 0; i < tree_paths.size(); ++i) {
          QueryVertexID path_root_vertex_idx = 0;
          for (QueryVertexID j = 0; j < tree_paths[i].size(); ++j) {
            VertexID cur_vertex = tree_paths[i][j];

            if (visited_vertices[cur_vertex]) continue;

            path_root_vertex_idx = j == 0 ? j : j - 1;
            break;
          }

          double cur_value = paths_embededdings_num[i][path_root_vertex_idx] /
                             (double)(candidates_sets_[tree_paths[i][path_root_vertex_idx]].size());
          if (cur_value < min_value) {
            min_value = cur_value;
            selected_path_index = i;
          }
        }

        for (QueryVertexID i = 0; i < tree_paths[selected_path_index].size(); ++i) {
          if (visited_vertices[tree_paths[selected_path_index][i]]) continue;

          order[selected_vertices_count] = tree_paths[selected_path_index][i];
          selected_vertices_count += 1;
          visited_vertices[tree_paths[selected_path_index][i]] = true;
        }

        tree_paths.erase(tree_paths.begin() + selected_path_index);
        paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
      }
    }

    // Order the leaves.
    while (!leaves.empty()) {
      double min_value = std::numeric_limits<double>::max();
      QueryVertexID selected_leaf_index = 0;

      for (QueryVertexID i = 0; i < leaves.size(); ++i) {
        QueryVertexID vertex = leaves[i];
        double cur_value = candidates_sets_[vertex].size();

        if (cur_value < min_value) {
          min_value = cur_value;
          selected_leaf_index = i;
        }
      }

      if (!visited_vertices[leaves[selected_leaf_index]]) {
        order[selected_vertices_count] = leaves[selected_leaf_index];
        selected_vertices_count += 1;
        visited_vertices[leaves[selected_leaf_index]] = true;
      }
      leaves.erase(leaves.begin() + selected_leaf_index);
    }
    return order;
  }
};

}  // namespace circinus

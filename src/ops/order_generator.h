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
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/bipartite_graph.h"
#include "graph/graph.h"
#include "graph/partitioned_graph.h"
#include "graph/query_graph.h"
#include "ops/logical/filter/cfl_filter.h"
#include "ops/logical/filter/daf_filter.h"
#include "ops/logical/filter/filter.h"
#include "ops/logical/filter/tso_filter.h"
#include "utils/query_utils.h"

namespace circinus {
class OrderGenerator {
 private:
  struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const {
      auto hash1 = std::hash<T1>{}(p.first);
      auto hash2 = std::hash<T2>{}(p.second);
      return hash1 ^ hash2;
    }
  };
  std::unordered_map<std::pair<QueryVertexID, QueryVertexID>, BipartiteGraph, hash_pair> bg_map_;

  const GraphBase* data_graph_;
  const QueryGraph* query_graph_;
  const std::vector<CandidateSetView>& candidates_;
  const std::vector<VertexID>& candidate_sizes_;
  const GraphMetadata& metadata_;

 public:
  OrderGenerator(const GraphBase* data_graph, const GraphMetadata& metadata, const QueryGraph* query_graph,
                 const std::vector<CandidateSetView>& candidates, const std::vector<VertexID>& candidate_sizes)
      : data_graph_(data_graph),
        query_graph_(query_graph),
        candidates_(candidates),
        candidate_sizes_(candidate_sizes),
        metadata_(metadata) {}

  const BipartiteGraph* getBipartiteGraph(QueryVertexID v1, QueryVertexID v2) {
    std::pair<QueryVertexID, QueryVertexID> p(v1, v2);
    auto bg = bg_map_.find(p);
    if (bg == bg_map_.end()) {
      BipartiteGraph newbg(v1, v2);
      newbg.populateGraph(data_graph_, candidates_);
      auto res = bg_map_.insert({p, std::move(newbg)});
      bg = res.first;
    }
    return &(bg->second);
  }

  std::vector<QueryVertexID> getOrder(OrderStrategy order_strategy) {
    switch (order_strategy) {
    case OrderStrategy::None:
    case OrderStrategy::CFL:
      return getCFLOrder();
    case OrderStrategy::DAF:
      return getDAFOrder();
    case OrderStrategy::TSO:
      return getTSOOrder();
    case OrderStrategy::GQL:
      return getGQLOrder();
    }
    return getCFLOrder();  // default
  }

  QueryVertexID selectGQLStartVertex() {
    QueryVertexID start_vertex = 0;
    QueryVertexID qg_v_cnt = query_graph_->getNumVertices();
    for (QueryVertexID i = 1; i < qg_v_cnt; ++i) {
      QueryVertexID cur_vertex = i;
      size_t size1 = candidate_sizes_[cur_vertex], size2 = candidate_sizes_[start_vertex];
      if (size1 < size2) {
        start_vertex = cur_vertex;
      } else if (size1 == size2 &&
                 query_graph_->getVertexOutDegree(cur_vertex) > query_graph_->getVertexOutDegree(start_vertex)) {
        start_vertex = cur_vertex;
      }
    }
    return start_vertex;
  }

  void updateValidVertices(QueryVertexID query_vertex, std::vector<bool>& visited, std::vector<bool>& adjacent) {
    visited[query_vertex] = true;
    auto[nbrs, cnt] = query_graph_->getOutNeighbors(query_vertex);
    for (uint32_t i = 0; i < cnt; ++i) {
      adjacent[nbrs[i]] = true;
    }
  }

  std::vector<QueryVertexID> getGQLOrder() {
    QueryVertexID qg_v_cnt = query_graph_->getNumVertices();
    std::vector<bool> visited_vertices(qg_v_cnt, false);
    std::vector<bool> adjacent_vertices(qg_v_cnt, false);
    std::vector<QueryVertexID> order(qg_v_cnt);

    QueryVertexID start_vertex = selectGQLStartVertex();
    order[0] = start_vertex;
    updateValidVertices(start_vertex, visited_vertices, adjacent_vertices);

    for (QueryVertexID i = 1; i < qg_v_cnt; ++i) {
      QueryVertexID next_vertex = 0;
      QueryVertexID min_value = metadata_.getNumVertices() + 1;
      for (QueryVertexID j = 0; j < qg_v_cnt; ++j) {
        QueryVertexID cur_vertex = j;
        if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
          size_t cnt = candidate_sizes_[cur_vertex];
          if (cnt < min_value) {
            min_value = cnt;
            next_vertex = cur_vertex;
          } else if (cnt == min_value &&
                     query_graph_->getVertexOutDegree(cur_vertex) > query_graph_->getVertexOutDegree(next_vertex)) {
            next_vertex = cur_vertex;
          }
        }
      }
      updateValidVertices(next_vertex, visited_vertices, adjacent_vertices);
      order[i] = next_vertex;
    }
    return order;
  }

  // FIXME(tatiana): need to support adaptive ordering strategy
  std::vector<QueryVertexID> getDAFOrder() {
    auto logical_filter = LogicalDAFFilter(metadata_, query_graph_, candidate_sizes_);
    return logical_filter.getBfsOrder();
  }

  void estimatePathEmbeddsingsNum(const std::vector<QueryVertexID>& path,
                                  std::vector<size_t>& estimated_embeddings_num) {
    CHECK_GT(path.size(), 1);
    std::vector<size_t> parent;
    std::vector<size_t> children;

    size_t begin = path.size() - 2, end = path.size() - 1;

    estimated_embeddings_num.resize(path.size() - 1);
    auto last_edge = getBipartiteGraph(path[begin], path[end]);
    children.resize(last_edge->getNumVertices());

    size_t sum = 0;
    for (auto& v : candidates_[path[begin]]) {
      uint32_t offset = last_edge->getOffset(v);
      children[offset] = last_edge->getVertexOutDegree(v);
      sum += children[offset];
    }

    estimated_embeddings_num[begin] = sum;

    for (uint32_t i = begin; i >= 1; --i) {
      begin = i - 1;
      end = i;
      auto edge = getBipartiteGraph(path[begin], path[end]);
      parent.resize(edge->getNumVertices());

      sum = 0;
      CHECK_EQ(path[end], last_edge->getSourceId());
      for (auto& v : candidates_[path[begin]]) {
        size_t local_sum = 0;
        auto vertex_view = edge->getOutNeighborsWithHint(v);
        for (auto nbr : vertex_view) {
          local_sum += children[last_edge->getOffset(nbr)];
        }
        parent[edge->getOffset(v)] = local_sum;
        sum += local_sum;
      }

      estimated_embeddings_num[begin] = sum;
      parent.swap(children);

      last_edge = edge;
    }
  }

  QueryVertexID generateNoneTreeEdgesCount(const std::vector<TreeNode>& tree_node,
                                           const std::vector<QueryVertexID>& path) {
    auto non_tree_edge_count = query_graph_->getVertexOutDegree(path[0]) - tree_node[path[0]].children_.size();
    for (size_t i = 1; i < path.size(); ++i) {
      auto vertex = path[i];
      non_tree_edge_count += query_graph_->getVertexOutDegree(vertex) - tree_node[vertex].children_.size() - 1;
    }

    return non_tree_edge_count;
  }

  void generateRootToLeafPaths(const std::vector<TreeNode>& tree_node, QueryVertexID cur_vertex,
                               std::vector<QueryVertexID>& cur_path, std::vector<std::vector<QueryVertexID>>& paths) {
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

  std::vector<QueryVertexID> getTSOOrder() {
    auto logical_filter = LogicalTSOFilter(metadata_, query_graph_, candidate_sizes_);
    const std::vector<TreeNode>& tree = logical_filter.getTree();
    const std::vector<QueryVertexID>& dfs_order = logical_filter.getDfsOrder();

    QueryVertexID qg_v_cnt = query_graph_->getNumVertices();
    std::vector<std::vector<QueryVertexID>> paths;
    paths.reserve(qg_v_cnt);

    std::vector<QueryVertexID> single_path;
    single_path.reserve(qg_v_cnt);

    generateRootToLeafPaths(tree, dfs_order[0], single_path, paths);
    std::vector<std::pair<double, std::vector<QueryVertexID>*>> path_orders;
    for (auto& path : paths) {
      std::vector<size_t> estimated_embeddings_num;
      QueryVertexID non_tree_edges_count = generateNoneTreeEdgesCount(tree, path);
      estimatePathEmbeddsingsNum(path, estimated_embeddings_num);
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

  void generateLeaves(std::vector<QueryVertexID>& leaves) {
    for (QueryVertexID i = 0; i < query_graph_->getNumVertices(); ++i) {
      QueryVertexID cur_vertex = i;
      if (query_graph_->getVertexOutDegree(cur_vertex) == 1) {
        leaves.push_back(cur_vertex);
      }
    }
  }

  void generateCorePaths(const std::vector<TreeNode>& tree_node, QueryVertexID cur_vertex,
                         std::vector<QueryVertexID>& cur_core_path, std::vector<std::vector<QueryVertexID>>& core_paths,
                         const TwoCoreSolver& tcs) {
    const TreeNode& node = tree_node[cur_vertex];
    cur_core_path.push_back(cur_vertex);

    bool is_core_leaf = true;
    for (const auto& child : node.children_) {
      if (tcs.isInCore(child)) {
        generateCorePaths(tree_node, child, cur_core_path, core_paths, tcs);
        is_core_leaf = false;
      }
    }

    if (is_core_leaf) {
      core_paths.emplace_back(cur_core_path);
    }
    cur_core_path.pop_back();
  }

  void generateTreePaths(const std::vector<TreeNode>& tree_node, QueryVertexID cur_vertex,
                         std::vector<QueryVertexID>& cur_tree_path,
                         std::vector<std::vector<QueryVertexID>>& tree_paths) {
    const TreeNode& node = tree_node[cur_vertex];
    cur_tree_path.push_back(cur_vertex);

    bool is_tree_leaf = true;
    for (auto child : node.children_) {
      if (query_graph_->getVertexOutDegree(child) > 1) {
        generateTreePaths(tree_node, child, cur_tree_path, tree_paths);
        is_tree_leaf = false;
      }
    }

    if (is_tree_leaf && cur_tree_path.size() > 1) {
      tree_paths.emplace_back(cur_tree_path);
    }
    cur_tree_path.pop_back();
  }

  std::vector<QueryVertexID> getCFLOrder() {
    auto logical_filter = LogicalCFLFilter(metadata_, query_graph_, candidate_sizes_);
    const std::vector<TreeNode>& tree = logical_filter.getTree();
    const std::vector<QueryVertexID>& bfs_order = logical_filter.getBfsOrder();

    QueryVertexID qg_v_cnt = query_graph_->getNumVertices();
    QueryVertexID root_vertex = bfs_order[0];
    std::vector<QueryVertexID> order(qg_v_cnt);
    std::vector<bool> visited_vertices(qg_v_cnt, false);

    std::vector<std::vector<QueryVertexID>> core_paths;
    std::vector<std::vector<std::vector<QueryVertexID>>> forests;
    std::vector<QueryVertexID> leaves;

    generateLeaves(leaves);

    const auto& tcs = logical_filter.getTwoCoreSolver();
    if (tcs.isInCore(root_vertex)) {
      std::vector<QueryVertexID> temp_core_path;
      generateCorePaths(tree, root_vertex, temp_core_path, core_paths, tcs);
      for (QueryVertexID i = 0; i < qg_v_cnt; ++i) {
        QueryVertexID cur_vertex = i;
        if (tcs.isInCore(cur_vertex)) {
          std::vector<std::vector<QueryVertexID>> temp_tree_paths;
          std::vector<QueryVertexID> temp_tree_path;
          // FIXME(tatiana): here the tree paths can contain core vertices?
          generateTreePaths(tree, cur_vertex, temp_tree_path, temp_tree_paths);
          if (!temp_tree_paths.empty()) {
            forests.emplace_back(temp_tree_paths);
          }
        }
      }
    } else {
      std::vector<std::vector<QueryVertexID>> temp_tree_paths;
      std::vector<QueryVertexID> temp_tree_path;
      generateTreePaths(tree, root_vertex, temp_tree_path, temp_tree_paths);
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
      paths_non_tree_edge_num.reserve(core_paths.size());
      for (auto& path : core_paths) {
        QueryVertexID non_tree_edge_num = generateNoneTreeEdgesCount(tree, path);
        paths_non_tree_edge_num.push_back(non_tree_edge_num + 1);

        std::vector<size_t> path_embeddings_num;
        estimatePathEmbeddsingsNum(path, path_embeddings_num);
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

      // Select subsequent paths
      // TODO(tatiana): the following codes and the codes for path ordering of forest can be shared?
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
                             (double)(candidate_sizes_[core_paths[i][path_root_vertex_idx]]);
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
    // TODO(tatiana): there is only path ordering within each tree, but not ordering trees in the forest?
    for (auto& tree_paths : forests) {
      std::vector<std::vector<size_t>> paths_embededdings_num;
      for (auto& path : tree_paths) {
        std::vector<size_t> path_embeddings_num;
        // FIXME(tatiana): no need to estimate for the path segment which contain query vertices in core? see line 335
        estimatePathEmbeddsingsNum(path, path_embeddings_num);
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

            // TODO(tatiana): j cannot be 0 because the first vertex of the path is in core? CHECK_NE(j, 0);
            path_root_vertex_idx = j == 0 ? j : j - 1;
            break;
          }

          double cur_value = paths_embededdings_num[i][path_root_vertex_idx] /
                             (double)(candidate_sizes_[tree_paths[i][path_root_vertex_idx]]);
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
    while (!leaves.empty()) {  // TODO(tatiana): simply sort the leaves instead of taking min and erasing?
      double min_value = std::numeric_limits<double>::max();
      QueryVertexID selected_leaf_index = 0;

      for (QueryVertexID i = 0; i < leaves.size(); ++i) {
        QueryVertexID vertex = leaves[i];
        double cur_value = candidate_sizes_[vertex];

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

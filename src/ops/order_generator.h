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
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "algorithms/partial_order.h"
#include "graph/bipartite_graph.h"
#include "graph/graph.h"
#include "graph/partitioned_graph.h"
#include "graph/query_graph.h"
#include "ops/logical/filter/cfl_filter.h"
#include "ops/logical/filter/daf_filter.h"
#include "ops/logical/filter/filter.h"
#include "ops/logical/filter/tso_filter.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class OrderGenerator {
 private:
  struct hash_pair {
    template <class T1, class T2>
    uint32_t operator()(const std::pair<T1, T2>& p) const {
      auto hash1 = std::hash<T1>{}(p.first);
      auto hash2 = std::hash<T2>{}(p.second);
      return hash1 ^ hash2;
    }
  };
  std::unordered_map<std::pair<QueryVertexID, QueryVertexID>, BipartiteGraph, hash_pair> bg_map_;

  const GraphBase* data_graph_ = nullptr;
  const QueryGraph* query_graph_;
  const std::vector<CandidateSetView>* candidates_ = nullptr;
  const std::vector<VertexID>* candidate_sizes_ = nullptr;
  const GraphMetadata* metadata_ = nullptr;

 public:
  OrderGenerator(const GraphBase* data_graph, const GraphMetadata& metadata, const QueryGraph* query_graph,
                 const std::vector<CandidateSetView>& candidates, const std::vector<VertexID>& candidate_sizes)
      : data_graph_(data_graph),
        query_graph_(query_graph),
        candidates_(&candidates),
        candidate_sizes_(&candidate_sizes),
        metadata_(&metadata) {}

  explicit OrderGenerator(const QueryGraph* query_graph) : query_graph_(query_graph) {}

  const BipartiteGraph* getBipartiteGraph(QueryVertexID v1, QueryVertexID v2) {
    std::pair<QueryVertexID, QueryVertexID> p(v1, v2);
    auto bg = bg_map_.find(p);
    if (bg == bg_map_.end()) {
      BipartiteGraph newbg(v1, v2);
      newbg.populateGraph(data_graph_, *candidates_);
      auto res = bg_map_.insert({p, std::move(newbg)});
      bg = res.first;
    }
    return &(bg->second);
  }

  std::vector<QueryVertexID> getOrder(OrderStrategy order_strategy, QueryVertexID seed_qv,
                                      const PartialOrder* po = nullptr) {
    switch (order_strategy) {
    case OrderStrategy::None:
      return getBFSOrder(seed_qv);
    case OrderStrategy::CFL:
      return getCFLOrder(seed_qv);
    case OrderStrategy::DAF:
      return getDAFOrder();
    case OrderStrategy::TSO:
      return getTSOOrder();
    case OrderStrategy::GQL:
      return getGQLOrder();
    case OrderStrategy::Online:
      return getOnlineOrder(seed_qv, po);
    }
    return getCFLOrder(seed_qv);  // default
  }

  QueryVertexID selectGQLStartVertex() {
    QueryVertexID start_vertex = 0;
    QueryVertexID qg_v_cnt = query_graph_->getNumVertices();
    for (QueryVertexID i = 1; i < qg_v_cnt; ++i) {
      QueryVertexID cur_vertex = i;
      uint32_t size1 = (*candidates_)[cur_vertex].size(), size2 = (*candidates_)[start_vertex].size();
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
      QueryVertexID min_value = (*metadata_).getNumVertices() + 1;
      for (QueryVertexID j = 0; j < qg_v_cnt; ++j) {
        QueryVertexID cur_vertex = j;
        if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
          uint32_t cnt = (*candidates_)[cur_vertex].size();
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
    auto logical_filter = LogicalDAFFilter((*metadata_), query_graph_, *candidate_sizes_);
    return logical_filter.getBfsOrder();
  }

  /** @returns True if no edge in any BipartiteGraph along the path */
  bool estimatePathEmbeddsingsNum(const std::vector<QueryVertexID>& path,
                                  std::vector<uint32_t>& estimated_embeddings_num) {
    CHECK_GT(path.size(), 1);
    std::vector<uint32_t> parent;
    std::vector<uint32_t> children;

    uint32_t begin = path.size() - 2, end = path.size() - 1;

    estimated_embeddings_num.resize(path.size() - 1);
    auto last_edge = getBipartiteGraph(path[begin], path[end]);
    if (last_edge->getNumEdges() == 0) {
      DLOG(INFO) << "no edge for BipartiteGraph between " << path[begin] << " and " << path[end] << " "
                 << (*candidates_)[begin].size() << " to " << (*candidates_)[end].size();
      return true;
    }
    children.resize(last_edge->getNumVertices());

    uint32_t sum = 0;
    for (auto& v : (*candidates_)[path[begin]]) {
      uint32_t offset = last_edge->getOffset(v);
      children[offset] = last_edge->getVertexOutDegree(v);
      sum += children[offset];
    }

    estimated_embeddings_num[begin] = sum;

    for (uint32_t i = begin; i >= 1; --i) {
      begin = i - 1;
      end = i;
      auto edge = getBipartiteGraph(path[begin], path[end]);
      if (edge->getNumEdges() == 0) {
        DLOG(INFO) << "no edge for BipartiteGraph between " << path[begin] << " and " << path[end] << " "
                   << (*candidates_)[begin].size() << " to " << (*candidates_)[end].size();
        return true;
      }
      parent.resize(edge->getNumVertices());

      sum = 0;
      CHECK_EQ(path[end], last_edge->getSourceId());
      for (auto& v : (*candidates_)[path[begin]]) {
        uint32_t local_sum = 0;
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
    return false;
  }

  QueryVertexID generateNoneTreeEdgesCount(const std::vector<TreeNode>& tree_node,
                                           const std::vector<QueryVertexID>& path) {
    auto non_tree_edge_count = query_graph_->getVertexOutDegree(path[0]) - tree_node[path[0]].children_.size();
    for (uint32_t i = 1; i < path.size(); ++i) {
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
    auto logical_filter = LogicalTSOFilter((*metadata_), query_graph_, *candidate_sizes_);
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
      std::vector<uint32_t> estimated_embeddings_num;
      QueryVertexID non_tree_edges_count = generateNoneTreeEdgesCount(tree, path);
      if (estimatePathEmbeddsingsNum(path, estimated_embeddings_num)) return {};
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
      if (verbosePlannerLog()) {
        LOG(INFO) << "[CFL] core path" << toString(cur_core_path);
      }
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
      if (verbosePlannerLog()) {
        LOG(INFO) << "[CFL] tree path" << toString(cur_tree_path);
      }
      tree_paths.emplace_back(cur_tree_path);
    }
    cur_tree_path.pop_back();
  }

  void SelectSubsequentPaths(std::vector<std::vector<QueryVertexID>>& paths,
                             std::vector<std::vector<uint32_t>>& paths_embededdings_num,
                             std::vector<bool>& visited_vertices, std::vector<QueryVertexID>& order,
                             uint32_t& selected_vertices_count) {
    while (!paths.empty()) {
      double min_value = std::numeric_limits<double>::max();
      uint32_t selected_path_index = 0;

      for (uint32_t i = 0; i < paths.size(); ++i) {
        QueryVertexID path_root_vertex_idx = 0;
        for (QueryVertexID j = 0; j < paths[i].size(); ++j) {
          QueryVertexID cur_vertex = paths[i][j];
          if (visited_vertices[cur_vertex]) continue;

          path_root_vertex_idx = j - 1;
          break;
        }

        if (paths_embededdings_num[i][path_root_vertex_idx] == 0) {
          LOG(WARNING) << "path score is 0, path starts from " << paths[i][path_root_vertex_idx];
        }
        double cur_value = paths_embededdings_num[i][path_root_vertex_idx] /
                           (double)((*candidates_)[paths[i][path_root_vertex_idx]].size());
        if (cur_value < min_value) {
          min_value = cur_value;
          selected_path_index = i;
        }
      }

      for (QueryVertexID i = 1; i < paths[selected_path_index].size(); ++i) {
        if (visited_vertices[paths[selected_path_index][i]]) continue;

        order[selected_vertices_count] = paths[selected_path_index][i];
        selected_vertices_count += 1;
        visited_vertices[paths[selected_path_index][i]] = true;
      }

      paths.erase(paths.begin() + selected_path_index);
      paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
    }
  }

  inline std::vector<QueryVertexID> getBFSOrder(QueryVertexID seed_qv) {
    auto logical_filter = LogicalCFLFilter(query_graph_, seed_qv);
    return logical_filter.getBfsOrder();
  }

  std::vector<QueryVertexID> getOnlineOrder(QueryVertexID seed_qv, const PartialOrder* po) {
    /* degree sorting to choose vertices with larger degree first */
    /* if degrees are equal, choose vertices with more constraints first */
    auto n_qvs = query_graph_->getNumVertices();
    std::vector<uint32_t> constraint_count;
    if (po != nullptr) {
      constraint_count = po->n_all_related_constraints;
    } else {
      constraint_count.resize(n_qvs, 0);
    }
    if (seed_qv == DUMMY_QUERY_VERTEX) {  // select start vertex
      seed_qv = 0;
      for (QueryVertexID v = 1; v < n_qvs; ++v) {
        if (query_graph_->getVertexOutDegree(v) > query_graph_->getVertexOutDegree(seed_qv) ||
            (query_graph_->getVertexOutDegree(v) == query_graph_->getVertexOutDegree(seed_qv) &&
             constraint_count[v] > constraint_count[seed_qv])) {
          seed_qv = v;
        }
      }
    }
    // degree and constraint sorting
    auto comp = [&constraint_count, q = query_graph_ ](QueryVertexID a, QueryVertexID b) {
      return q->getVertexOutDegree(a) < q->getVertexOutDegree(b) ||
             (q->getVertexOutDegree(a) == q->getVertexOutDegree(b) && constraint_count[a] < constraint_count[b]);
    };
    std::vector<QueryVertexID> matching_order;
    std::priority_queue<QueryVertexID, std::vector<QueryVertexID>, decltype(comp)> queue(comp);
    std::vector<bool> visited(n_qvs, 0);
    matching_order.reserve(n_qvs);
    queue.push(seed_qv);
    visited[seed_qv] = true;
    while (!queue.empty()) {
      auto v = queue.top();
      queue.pop();
      matching_order.push_back(v);
      auto nbrs = query_graph_->getOutNeighbors(v);
      for (uint32_t i = 0; i < nbrs.second; ++i) {
        if (visited[nbrs.first[i]]) continue;
        queue.push(nbrs.first[i]);
        visited[nbrs.first[i]] = true;
      }
    }
    CHECK_EQ(matching_order.size(), query_graph_->getNumVertices());
    return matching_order;
  }

  std::vector<QueryVertexID> getCFLOrder(QueryVertexID seed_qv) {
    auto logical_filter = LogicalCFLFilter((*metadata_), query_graph_, *candidate_sizes_, seed_qv);
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
    if (verbosePlannerLog()) {
      LOG(INFO) << "[CFL] leaves" << toString(leaves);
    }

    const auto& tcs = logical_filter.getTwoCoreSolver();
    if (tcs.isInCore(root_vertex)) {
      std::vector<QueryVertexID> temp_core_path;
      generateCorePaths(tree, root_vertex, temp_core_path, core_paths, tcs);
      for (QueryVertexID i = 0; i < qg_v_cnt; ++i) {
        QueryVertexID cur_vertex = i;
        if (tcs.isInCore(cur_vertex)) {
          std::vector<std::vector<QueryVertexID>> temp_tree_paths;
          std::vector<QueryVertexID> temp_tree_path;
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
      std::vector<std::vector<uint32_t>> paths_embededdings_num;
      std::vector<QueryVertexID> paths_non_tree_edge_num;
      paths_non_tree_edge_num.reserve(core_paths.size());
      for (auto& path : core_paths) {
        QueryVertexID non_tree_edge_num = generateNoneTreeEdgesCount(tree, path);
        paths_non_tree_edge_num.push_back(non_tree_edge_num + 1);

        std::vector<uint32_t> path_embeddings_num;
        if (estimatePathEmbeddsingsNum(path, path_embeddings_num)) return {};
        paths_embededdings_num.emplace_back(std::move(path_embeddings_num));
      }

      // Select the start path.
      double min_value = std::numeric_limits<double>::max();
      uint32_t selected_path_index = 0;

      for (uint32_t i = 0; i < core_paths.size(); ++i) {
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

      SelectSubsequentPaths(core_paths, paths_embededdings_num, visited_vertices, order, selected_vertices_count);
    }

    // Order tree paths.
    // TODO(tatiana): there is only path ordering within each tree, but not ordering trees in the forest?
    for (auto& tree_paths : forests) {
      std::vector<std::vector<uint32_t>> paths_embededdings_num;
      for (auto& path : tree_paths) {
        std::vector<uint32_t> path_embeddings_num;
        // FIXME(tatiana): no need to estimate for the path segment which contain query vertices in core?
        if (estimatePathEmbeddsingsNum(path, path_embeddings_num)) return {};
        paths_embededdings_num.emplace_back(path_embeddings_num);
      }
      SelectSubsequentPaths(tree_paths, paths_embededdings_num, visited_vertices, order, selected_vertices_count);
    }

    // Order the leaves.
    sort(leaves.begin(), leaves.end(),
         [&](QueryVertexID l, QueryVertexID r) { return (*candidates_)[l].size() < (*candidates_)[r].size(); });
    for (QueryVertexID leaf : leaves) {
      if (!visited_vertices[leaf]) {
        order[selected_vertices_count++] = leaf;  // no need to update visited_vertices
      }
    }
    CHECK_EQ(selected_vertices_count, qg_v_cnt);
    return order;
  }
};

}  // namespace circinus

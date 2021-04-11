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

#include "plan/naive_planner.h"

#include <algorithm>
#include <bitset>
#include <numeric>
#include <queue>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "algorithms/k_core.h"
#include "algorithms/minimum_weight_vertex_cover.h"
#include "algorithms/vertex_cover.h"
#include "graph/query_graph.h"
#include "plan/execution_plan.h"
#include "utils/hashmap.h"

namespace circinus {

bool NaivePlanner::hasValidCandidate() {
  uint32_t i = 0;
  for (auto& cardinality : *candidate_cardinality_) {
    if (cardinality < 0) {
      LOG(INFO) << "candidate for vertex is empty " << i;
      return false;
    }
    ++i;
  }
  return true;
}

ExecutionPlan* NaivePlanner::generatePlan(const std::vector<QueryVertexID>& use_order, Profiler* profiler) {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return nullptr;
  }
  if (use_order.empty()) {
    auto log_cardinality = logCardinality();
    WeightedBnB vertex_cover_solver(query_graph_, log_cardinality);
    vertex_cover_solver.computeVertexCover();
    auto& covers = vertex_cover_solver.getBestCovers();
    CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
    auto select_cover = covers[BnB::getSmallestCover(covers).first.front()];
    std::stringstream ss;
    std::string s = "";
    for (QueryVertexID i = 0; i < select_cover.size(); ++i) {
      ss << ' ' << i << ':' << (select_cover[i] == 1 ? "key" : "set");
      s += (select_cover[i] == 1) ? "1," : "0,";
    }
    LOG(INFO) << "cover" << ss.str();

    double res = 1;
    for (uint32_t i = 0; i < select_cover.size(); ++i) {
      if (select_cover[i] == 1) {
        res *= (*candidate_cardinality_)[i];
      }
    }
    LOG(INFO) << res << "  " << s;

    TwoCoreSolver solver;
    auto& core_table = solver.get2CoreTable(query_graph_);
    // now we only consider a random smallest MWVC from covers
    QueryVertexID v = 0;
    std::vector<QueryVertexID> cover;
    for (auto assignment : select_cover) {
      if (assignment == 1) {
        cover.push_back(v);
      }
      ++v;
    }
    // start matching from the vertex with the smallest cardinality in cover
    auto start_vertex = selectStartingVertex(cover);

    matching_order_ = generateMatchingOrder(query_graph_, core_table, start_vertex);
  } else {
    matching_order_ = use_order;
  }
  plan_.populatePhysicalPlan(query_graph_, matching_order_, getCoverByOrder(), profiler);
  return &plan_;
}

ExecutionPlan* NaivePlanner::generatePlanWithoutCompression(const std::vector<QueryVertexID>& use_order,
                                                            Profiler* profiler) {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return nullptr;
  }
  // all query vertices are in cover, so that no compression is done
  std::vector<int> select_cover(query_graph_->getNumVertices(), 1);

  if (use_order.empty()) {
    TwoCoreSolver solver;
    auto& core_table = solver.get2CoreTable(query_graph_);
    // now we only consider a random smallest MWVC from covers
    std::vector<QueryVertexID> cover(query_graph_->getNumVertices());
    std::iota(cover.begin(), cover.end(), 0);
    // start matching from the vertex with the smallest cardinality in cover
    auto start_vertex = selectStartingVertex(cover);
    matching_order_ = generateMatchingOrder(query_graph_, core_table, start_vertex);
  } else {
    matching_order_ = use_order;
  }

  plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover, profiler);
  return &plan_;
}

std::vector<QueryVertexID> NaivePlanner::generateMatchingOrder(const QueryGraph* g, const std::vector<int>& core_table,
                                                               QueryVertexID start_vertex) {
  std::vector<bool> visited(g->getNumVertices(), false);  // visited when pushed into the queue
  std::vector<QueryVertexID> order(g->getNumVertices());  // the order out of queue
  // if both in core or both not , prioritize the vertex with smaller cardinality; otherwise prioritize the core
  auto compare = [this, &core_table](QueryVertexID v1, QueryVertexID v2) {
    if (TwoCoreSolver::isInCore(core_table, v1)) {
      if (TwoCoreSolver::isInCore(core_table, v2)) {
        return (*candidate_cardinality_)[v1] > (*candidate_cardinality_)[v2];
      }
      return false;
    }
    if (TwoCoreSolver::isInCore(core_table, v2)) {
      return true;
    }
    return (*candidate_cardinality_)[v1] > (*candidate_cardinality_)[v2];
  };
  std::priority_queue<QueryVertexID, std::vector<QueryVertexID>, decltype(compare)> queue(compare);
  // initialization
  visited[start_vertex] = true;
  order.front() = start_vertex;
  uint32_t order_offset = 1;
  auto neighbors = g->getOutNeighbors(start_vertex);
  for (uint32_t i = 0; i < neighbors.second; ++i) {
    visited[neighbors.first[i]] = true;
    queue.push(neighbors.first[i]);
  }
  // repeatedly traverse to a best vertex
  while (!queue.empty()) {
    auto next = queue.top();
    queue.pop();
    order[order_offset++] = next;
    auto neighbors = g->getOutNeighbors(next);
    for (uint32_t i = 0; i < neighbors.second; ++i) {
      if (!visited[neighbors.first[i]]) {
        visited[neighbors.first[i]] = true;
        queue.push(neighbors.first[i]);
      }
    }
  }
  return order;
}

QueryVertexID NaivePlanner::selectStartingVertex(const std::vector<QueryVertexID>& cover) {
  // std::sort(cover.begin(), cover.end(), [this](QueryVertexID v1, QueryVertexID v2) { return
  // (*candidate_cardinality_)[v1] / query_graph_->getVertexOutDegree(v1) < (*candidate_cardinality_)[v2] /
  // query_graph_->getVertexOutDegree(v2); });

  QueryVertexID start_vertex = cover.front();
  double min = (*candidate_cardinality_)[cover.front()];

  for (auto v : cover) {
    double score = (*candidate_cardinality_)[v] / query_graph_->getVertexOutDegree(v);
    DLOG(INFO) << "cover vertex " << v << " score " << score << " min " << min;
    if (score < min) {
      min = score;
      start_vertex = v;
    }
  }
  return start_vertex;
}

unordered_map<QueryVertexID, uint32_t> NaivePlanner::getDynamicCoreCoverEager(const std::vector<int>& select_cover) {
  unordered_map<QueryVertexID, uint32_t> level_become_key;
  DCHECK(hasValidCandidate());
  auto cover = select_cover;
  unordered_set<QueryVertexID> existing_vertices;
  std::vector<uint32_t> vertex_order(matching_order_.size());
  existing_vertices.reserve(matching_order_.size());
  existing_vertices.insert(matching_order_.front());
  vertex_order[matching_order_.front()] = 0;
  // for the first query vertex, we always put it as key if it is in the cover. if this is changed, the implementation
  // for populating execution plan must be changed accordingly
  if (cover[matching_order_.front()] == 1) {
    level_become_key.insert({matching_order_.front(), 0});
  }
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    auto v = matching_order_[i];
    vertex_order[v] = i;
    auto nb = query_graph_->getOutNeighbors(v);
    if (cover[v] == 1) {  // vertex v is in key, we eagerly check whether the current target can be put in key later
      bool delay_enumeration = true;
      // can be put in key later for a larger subquery if all the existing vertices connecting to v is in key
      for (uint32_t j = 0; j < nb.second; ++j) {
        if (existing_vertices.count(nb.first[j]) && cover[nb.first[j]] != 1) {
          delay_enumeration = false;
          break;
        }
      }
      if (!delay_enumeration) {
        level_become_key.insert({v, i});  // the level is the subquery size - 1
        cover[v] = 1;
      } else {
        cover[v] = 0;
      }
    } else {  // vertex v is in set, we need to ensure that all its neighbors are in key
      for (uint32_t j = 0; j < nb.second; ++j) {  // restore the delayed vertex as key
        if (existing_vertices.count(nb.first[j]) && cover[nb.first[j]] != 1) {
          CHECK_EQ(select_cover[nb.first[j]], 1);
          cover[nb.first[j]] = 1;
          if ((i - vertex_order[nb.first[j]]) == 1) {  // if enumerate in the consecutive subquery, then do not delay
            level_become_key.insert({nb.first[j], vertex_order[nb.first[j]]});
          } else {
            level_become_key.insert({nb.first[j], i});  // the level is the subquery size - 1
          }
        }
      }
    }
    existing_vertices.insert(v);
  }
  return level_become_key;
}

std::pair<uint32_t, uint32_t> NaivePlanner::analyzeDynamicCoreCoverEagerInner(const std::vector<int>& select_cover) {
  DCHECK(hasValidCandidate());
  uint32_t sum_delayed_steps = 0;
  auto cover = select_cover;
  unordered_set<QueryVertexID> existing_vertices;
  std::vector<uint32_t> vertex_order(matching_order_.size());
  existing_vertices.reserve(matching_order_.size());
  existing_vertices.insert(matching_order_.front());
  vertex_order[matching_order_.front()] = 0;
  auto current_subquery_key_size = (select_cover[matching_order_.front()] == 1);
  uint32_t key_sizes = current_subquery_key_size;
  LOG(INFO) << "analyzeDynamicCoreCoverEagerInner";
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    auto v = matching_order_[i];
    vertex_order[v] = i;
    auto nb = query_graph_->getOutNeighbors(v);
    if (cover[v] == 1) {  // vertex v is in key, we eagerly check whether the current target can be put in key later
      bool delay_enumeration = true;
      // can be put in key later for a larger subquery if all the existing vertices connecting to v is in key
      for (uint32_t j = 0; j < nb.second; ++j) {
        if (existing_vertices.count(nb.first[j]) && cover[nb.first[j]] != 1) {
          delay_enumeration = false;
          break;
        }
      }
      cover[v] = !delay_enumeration;
      current_subquery_key_size += !delay_enumeration;
    } else {  // vertex v is in set, we need to ensure that all its neighbors are in key
      for (uint32_t j = 0; j < nb.second; ++j) {  // restore the delayed vertex as key
        if (existing_vertices.count(nb.first[j]) && cover[nb.first[j]] != 1) {
          CHECK_EQ(select_cover[nb.first[j]], 1);
          cover[nb.first[j]] = 1;
          current_subquery_key_size += 1;
          if ((i - vertex_order[nb.first[j]]) == 1) {  // if enumerate in the consecutive subquery, then do not delay
            current_subquery_key_size += 1;            // restore for the last subquery
          } else {
            sum_delayed_steps += (i - vertex_order[nb.first[j]]);
            LOG(INFO) << nb.first[j] << " delayed from " << vertex_order[nb.first[j]] << " for "
                      << (i - vertex_order[nb.first[j]]) << " steps";
          }
        }
      }
    }
    key_sizes += current_subquery_key_size;
    existing_vertices.insert(v);
  }
  return std::make_pair(key_sizes, sum_delayed_steps);
}

std::pair<uint32_t, uint32_t> NaivePlanner::analyzeDynamicCoreCoverEager(const std::vector<QueryVertexID>& use_order) {
  CHECK(hasValidCandidate());
  auto log_cardinality = logCardinality();
  WeightedBnB vertex_cover_solver(query_graph_, log_cardinality);
  vertex_cover_solver.computeVertexCover();
  auto& covers = vertex_cover_solver.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto select_cover = covers[BnB::getSmallestCover(covers).first.front()];
  if (!use_order.empty()) {
    matching_order_ = use_order;
  }
  if (matching_order_.empty()) {
    TwoCoreSolver solver;
    auto& core_table = solver.get2CoreTable(query_graph_);
    // now we only consider a random smallest MWVC from covers
    QueryVertexID v = 0;
    std::vector<QueryVertexID> cover;
    for (auto assignment : select_cover) {
      if (assignment == 1) {
        cover.push_back(v);
      }
      ++v;
    }
    // start matching from the vertex with the smallest cardinality in cover
    auto start_vertex = selectStartingVertex(cover);
    matching_order_ = generateMatchingOrder(query_graph_, core_table, start_vertex);
  }
  return analyzeDynamicCoreCoverEagerInner(select_cover);
}

ExecutionPlan* NaivePlanner::generatePlanWithEagerDynamicCover(const std::vector<QueryVertexID>& use_order,
                                                               Profiler* profiler) {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return nullptr;
  }
  auto log_cardinality = logCardinality();
  WeightedBnB vertex_cover_solver(query_graph_, log_cardinality);
  vertex_cover_solver.computeVertexCover();
  auto& covers = vertex_cover_solver.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto smallest_covers = BnB::getSmallestCover(covers);
  auto select_cover = covers[smallest_covers.first.front()];

  if (use_order.empty()) {
    TwoCoreSolver solver;
    auto& core_table = solver.get2CoreTable(query_graph_);
    // now we only consider a random smallest MWVC from covers
    QueryVertexID v = 0;
    std::vector<QueryVertexID> cover;
    for (auto assignment : select_cover) {
      if (assignment == 1) {
        cover.push_back(v);
      }
      ++v;
    }
    // start matching from the vertex with the smallest cardinality in cover
    auto start_vertex = selectStartingVertex(cover);
    matching_order_ = generateMatchingOrder(query_graph_, core_table, start_vertex);
  } else {
    matching_order_ = use_order;
  }

  // select a cover with minimum key sizes
  auto select_cover_index = 0;
  uint32_t min_key_sizes = ~0u;
  for (auto cover_idx : smallest_covers.first) {
    auto current_key_sizes = analyzeDynamicCoreCoverEagerInner(covers[cover_idx]).first;
    if (current_key_sizes < min_key_sizes) {
      min_key_sizes = current_key_sizes;
      select_cover_index = cover_idx;
    }
  }

  plan_.populatePhysicalPlan(query_graph_, matching_order_, covers[select_cover_index],
                             getDynamicCoreCoverEager(covers[select_cover_index]));
  plan_.setProfiler(profiler);
  return &plan_;
}

ExecutionPlan* NaivePlanner::generatePlanWithDynamicCover(Profiler* profiler) {
  std::vector<std::vector<double>> costs_car(covers_.size());
  std::vector<std::vector<uint32_t>> pre(covers_.size());
  std::vector<std::vector<double>> car(covers_.size());

  costs_car[0].resize(covers_[0].size(), 0);
  for (uint32_t i = 0; i < covers_.size(); ++i) {
    car[i].resize(covers_[i].size());
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      car[i][j] = 1;
      for (QueryVertexID qid : covers_[i][j].cover) {
        car[i][j] = car[i][j] * (*candidate_cardinality_)[qid];
      }
    }
  }
  for (uint32_t i = 1; i < covers_.size(); ++i) {
    costs_car[i].resize(covers_[i].size(), -1);
    pre[i].resize(covers_[i].size());
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      auto& cover_node = covers_[i][j];
      for (uint32_t par : cover_node.parents) {
        if (costs_car[i - 1][par] >= 0) {
          if (costs_car[i][j] < 0 || costs_car[i - 1][par] + car[i - 1][par] < costs_car[i][j]) {
            costs_car[i][j] = costs_car[i - 1][par] + car[i - 1][par];
            pre[i][j] = par;
          }
        }
      }
    }
  }

  CHECK_GT(covers_.size(), 0);
  uint32_t last = covers_.size() - 1;
  double mini_cost = -1;
  int best_idx = -1;
  for (uint32_t i = 0; i < covers_[last].size(); ++i) {
    std::vector<int> select_cover(matching_order_.size(), 0);
    for (uint32_t j = 0; j < matching_order_.size(); ++j) {
      if (covers_[last][i].cover_bits >> j & 1) {
        select_cover[j] = 1;
      }
    }

    std::string s = "";
    for (uint32_t j = 0; j < matching_order_.size(); ++j) {
      s += std::to_string(select_cover[j]) + ",";
    }
    LOG(INFO) << i << " cover " << s << " cardinality " << car[last][i] << ", costs_car " << costs_car[last][i];

    if (costs_car[last][i] >= 0 && (mini_cost < 0 || mini_cost > costs_car[last][i])) {
      mini_cost = costs_car[last][i];
      best_idx = i;
    }
  }
  CHECK(best_idx != -1) << "ERROR: can not get best plan index.";
  // best_idx = 0;
  LOG(INFO) << "best last level cover idx " << best_idx << "  " << car[last][best_idx];
  std::vector<int> select_cover(matching_order_.size(), 0);
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    if (covers_[last][best_idx].cover_bits >> i & 1) {
      select_cover[i] = 1;
    }
  }

  LOG(INFO) << mini_cost;
  std::string s = "";
  for (uint32_t i = 0; i < select_cover.size(); ++i) {
    s += std::to_string(select_cover[i]) + " ";
  }
  DLOG(INFO) << s;

  unordered_map<QueryVertexID, uint32_t> level_become_key;
  std::vector<uint32_t> best_path;
  best_path.emplace_back(best_idx);
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    std::string s = "[";
    for (auto vid : covers_[last][best_idx].cover) {
      s += std::to_string(vid) + " ";
    }
    s += "]";
    DLOG(INFO) << s << " " << covers_[last][best_idx].cover_bits;
    best_idx = pre[last][best_idx];
    best_path.emplace_back(best_idx);
    last--;
  }
  s = "[";
  for (auto vid : covers_[last][best_idx].cover) {
    s += std::to_string(vid) + " ";
  }
  s += "]";
  DLOG(INFO) << s << " " << covers_[last][best_idx].cover_bits;
  std::reverse(best_path.begin(), best_path.end());

  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    for (auto vid : covers_[i][best_path[i]].cover) {
      if (vid != matching_order_[i] && i > 0 && !(covers_[i - 1][best_path[i - 1]].cover_bits >> vid & 1)) {
        // if the first vertex becomes key in the next subquery, make it as key initially
        level_become_key.insert({vid, (i == 1) ? 0 : i});
        DLOG(INFO) << vid << " " << i;
      }
      if (vid == matching_order_[i]) {
        level_become_key.insert({vid, i});
        DLOG(INFO) << vid << " " << i;
      }
    }
  }
  plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover, level_become_key);
  plan_.setProfiler(profiler);

  return &plan_;
}

ExecutionPlan* NaivePlanner::generatePlanWithSampleExecution(const std::vector<std::vector<double>>& cardinality,
                                                             const std::vector<double>& level_cost,
                                                             Profiler* profiler) {
  std::vector<std::vector<double>> costs_car_sample_cost;
  std::vector<std::vector<double>> costs_car;
  std::vector<std::vector<uint32_t>> pre;
  costs_car_sample_cost.resize(covers_.size());
  costs_car.resize(covers_.size());
  pre.resize(covers_.size());
  costs_car_sample_cost[0].resize(covers_[0].size(), 0);
  costs_car[0].resize(covers_[0].size(), 0);
  std::vector<std::vector<double>> car;
  car.resize(covers_.size());
  for (uint32_t i = 0; i < covers_.size(); ++i) {
    car[i].resize(covers_[i].size());
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      car[i][j] = 1;
      for (QueryVertexID qid : covers_[i][j].cover) {
        car[i][j] = car[i][j] * cardinality[i][qid];
      }
    }
  }
  for (uint32_t i = 1; i < covers_.size(); ++i) {
    costs_car[i].resize(covers_[i].size(), -1);
    costs_car_sample_cost[i].resize(covers_[i].size(), -1);
    pre[i].resize(covers_[i].size());
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      auto& cover_node = covers_[i][j];
      for (uint32_t par : cover_node.parents) {
        if (costs_car_sample_cost[i - 1][par] >= 0) {
          if (costs_car_sample_cost[i][j] < 0 ||
              costs_car_sample_cost[i - 1][par] + car[i - 1][par] * level_cost[i - 1] < costs_car_sample_cost[i][j]) {
            costs_car_sample_cost[i][j] = costs_car_sample_cost[i - 1][par] + car[i - 1][par] * level_cost[i - 1];
            pre[i][j] = par;
          }
        }
        if (costs_car[i - 1][par] >= 0) {
          if (costs_car[i][j] < 0 || costs_car[i - 1][par] + car[i - 1][par] < costs_car[i][j]) {
            costs_car[i][j] = costs_car[i - 1][par] + car[i - 1][par];
            pre[i][j] = par;
          }
        }
      }
    }
  }

  CHECK_GT(covers_.size(), 0);
  uint32_t last = covers_.size() - 1;
  double mini_cost = -1;
  int best_idx = -1;
  for (uint32_t i = 0; i < covers_[last].size(); ++i) {
    std::vector<int> select_cover(matching_order_.size(), 0);
    for (uint32_t j = 0; j < matching_order_.size(); ++j) {
      select_cover[j] = covers_[last][i].cover_bits >> j & 1;
    }

    std::string s = "";
    for (uint32_t j = 0; j < matching_order_.size(); ++j) {
      s += std::to_string(select_cover[j]) + ",";
    }
    LOG(INFO) << i << " cover " << s << " cardinality " << car[last][i] << ", costs_car_sample_cost "
              << costs_car_sample_cost[last][i] << ", costs_car " << costs_car[last][i];

    if (costs_car_sample_cost[last][i] >= 0 && (mini_cost < 0 || mini_cost > costs_car_sample_cost[last][i])) {
      mini_cost = costs_car_sample_cost[last][i];
      best_idx = i;
    }
  }
  CHECK(best_idx != -1) << "ERROR: can not get best plan index.";
  // best_idx = 2;
  LOG(INFO) << "best last level cover idx " << best_idx << "  " << car[last][best_idx];
  std::vector<int> select_cover(matching_order_.size(), 0);
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    select_cover[i] = covers_[last][best_idx].cover_bits >> i & 1;
  }

  LOG(INFO) << mini_cost;
  std::string s = "";
  for (uint32_t i = 0; i < select_cover.size(); ++i) {
    s += std::to_string(select_cover[i]) + " ";
  }
  DLOG(INFO) << s;

  unordered_map<QueryVertexID, uint32_t> level_become_key;
  std::vector<uint32_t> best_path;
  best_path.emplace_back(best_idx);
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    std::string s = "[";
    for (auto vid : covers_[last][best_idx].cover) {
      s += std::to_string(vid) + " ";
    }
    s += "]";
    DLOG(INFO) << s << " " << covers_[last][best_idx].cover_bits;
    best_idx = pre[last][best_idx];
    best_path.emplace_back(best_idx);
    last--;
  }

  std::reverse(best_path.begin(), best_path.end());
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    for (auto vid : covers_[i][best_path[i]].cover) {
      if (vid != matching_order_[i] && !(covers_[i - 1][best_path[i - 1]].cover_bits >> vid & 1)) {
        if (i == 1) {
          level_become_key.insert({vid, i - 1});
          DLOG(INFO) << vid << " " << i - 1;
        } else {
          level_become_key.insert({vid, i});
          DLOG(INFO) << vid << " " << i;
        }
      }
      if (vid == matching_order_[i]) {
        level_become_key.insert({vid, i});
        DLOG(INFO) << vid << " " << i;
      }
    }
  }
  plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover, level_become_key);
  plan_.setProfiler(profiler);

  return &plan_;
}

void NaivePlanner::generateOrder(const std::vector<QueryVertexID>& use_order) {
  if (!hasValidCandidate()) {
    return;
  }

  // this is for getting matching order
  auto log_cardinality = logCardinality();
  WeightedBnB vertex_cover_solver(query_graph_, log_cardinality);
  vertex_cover_solver.computeVertexCover();
  auto& covers = vertex_cover_solver.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto smallest_covers = BnB::getSmallestCover(covers);
  auto select_cover = covers[smallest_covers.first.front()];
  std::string s = "";
  for (uint32_t j = 0; j < select_cover.size(); ++j) {
    s += std::to_string(select_cover[j]) + " ";
  }
  LOG(INFO) << s;

  if (use_order.empty()) {
    TwoCoreSolver solver;
    auto& core_table = solver.get2CoreTable(query_graph_);
    // now we only consider a random smallest MWVC from covers
    QueryVertexID v = 0;
    std::vector<QueryVertexID> cover;
    for (auto assignment : select_cover) {
      if (assignment == 1) {
        cover.push_back(v);
      }
      ++v;
    }
    // start matching from the vertex with the smallest cardinality in cover
    auto start_vertex = selectStartingVertex(cover);
    matching_order_ = generateMatchingOrder(query_graph_, core_table, start_vertex);
  } else {
    matching_order_ = use_order;
  }
}

std::vector<int> NaivePlanner::getCoverByOrder() const {
  // get a best cover for the query graph without the last matched vertex, the last matched vertex does not affect the
  // amount of intersection
  auto subquery_vertices = matching_order_;
  subquery_vertices.pop_back();
  auto subquery = query_graph_->getInducedSubgraph(subquery_vertices);

  std::vector<double> candidate_cardinality_by_match_order(matching_order_.size() - 1);
  for (uint32_t i = 0; i < candidate_cardinality_by_match_order.size(); ++i) {
    candidate_cardinality_by_match_order[i] = log2((*candidate_cardinality_)[matching_order_[i]]);
  }

  WeightedBnB vc_solver(&subquery, candidate_cardinality_by_match_order);
  vc_solver.computeVertexCover();
  auto& covers = vc_solver.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto& select_cover = covers[BnB::getSmallestCover(covers).first.front()];
  std::vector<int> ret(matching_order_.size());
  for (uint32_t i = 0; i < select_cover.size(); ++i) {
    ret[matching_order_[i]] = select_cover[i];
  }
  bool all_key = true;
  auto nbs = query_graph_->getOutNeighbors(matching_order_.back());
  for (uint32_t i = 0; i < nbs.second; ++i) {
    if (ret[nbs.first[i]] != 1) {
      all_key = false;
      break;
    }
  }
  ret[matching_order_.back()] = !all_key;
  std::stringstream ss;
  for (QueryVertexID i = 0; i < ret.size(); ++i) {
    ss << ' ' << i << ':' << (ret[i] == 1 ? "key" : "set");
  }
  LOG(INFO) << "cover" << ss.str();
  return ret;
}

void NaivePlanner::generateCoverNode(const std::vector<std::vector<double>>& cardinality) {
  if (!hasValidCandidate()) {
    return;
  }

  auto subquery_vertices = matching_order_;
  covers_.resize(matching_order_.size());

  std::vector<unordered_set<QueryVertexID>> existing_vertices;
  existing_vertices.resize(matching_order_.size());
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    for (uint32_t j = i + 1; j < matching_order_.size(); ++j) {
      existing_vertices[j].insert(matching_order_[i]);
    }
  }

  std::vector<std::vector<QueryVertexID>> to_intersect_vertices;
  to_intersect_vertices.resize(matching_order_.size());
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    QueryVertexID qid = matching_order_[i];
    auto nbrs = query_graph_->getOutNeighbors(qid);
    for (uint32_t j = 0; j < nbrs.second; ++j) {
      if (existing_vertices[i].find(nbrs.first[j]) != existing_vertices[i].end()) {
        to_intersect_vertices[i].emplace_back(nbrs.first[j]);
      }
    }
    std::sort(to_intersect_vertices[i].begin(), to_intersect_vertices[i].end(),
              [&](QueryVertexID qid1, QueryVertexID qid2) { return cardinality[i][qid1] > cardinality[i][qid2]; });
  }

  std::vector<std::vector<double>> candidate_cardinality_by_match_order(matching_order_.size());
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    candidate_cardinality_by_match_order[i].resize(i + 1);
    for (uint32_t j = 0; j < i + 1; ++j) {
      candidate_cardinality_by_match_order[i][j] = log2(cardinality[i][matching_order_[j]]);
    }
  }

  for (uint32_t i = matching_order_.size() - 1; i > 0; --i) {
    subquery_vertices.resize(i + 1);
    // compute a vertex cover of subquery i
    auto subquery = query_graph_->getInducedSubgraph(subquery_vertices);
    WeightedBnB vc_solver(&subquery, candidate_cardinality_by_match_order[i]);
    vc_solver.computeVertexCover();
    auto& covers = vc_solver.getBestCovers();
    CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
    auto choice = BnB::getSmallestCover(covers);
    auto select_cover = covers[choice.first.front()];
    CoverNode new_cover_node;
    new_cover_node.cover_bits = 0;
    for (uint32_t j = 0; j < select_cover.size(); ++j) {
      if (select_cover[j] == 1) {
        QueryVertexID v = matching_order_[j];
        new_cover_node.cover_bits |= 1ULL << v;
        new_cover_node.cover.push_back(v);
      }
    }

    if (i == matching_order_.size() - 1) {
      std::string s = "";
      std::vector<uint32_t> mapped_cover(select_cover.size(), 0);
      for (uint32_t j = 0; j < select_cover.size(); ++j) {
        if (select_cover[j] == 1) {
          mapped_cover[matching_order_[j]] = 1;
        }
      }
      for (uint32_t j = 0; j < mapped_cover.size(); ++j) {
        s += std::to_string(mapped_cover[j]) + ",";
      }
      DLOG(INFO) << "subquery " << i << " cover " << s;
    }

    covers_[i].emplace_back(new_cover_node);
    CoverNode nxt_cover_node = new_cover_node;
    for (uint32_t j = i + 1; j < matching_order_.size(); ++j) {
      QueryVertexID new_v = matching_order_[j];
      uint32_t to_intersect_vertices_all_key = 1;
      double set_cardinality = 1;
      for (QueryVertexID existing_v : to_intersect_vertices[j]) {
        if (!(nxt_cover_node.cover_bits >> existing_v & 1)) {
          to_intersect_vertices_all_key = 0;
          set_cardinality *= cardinality[j][existing_v];
        }
      }

      // we need to add key
      if (!to_intersect_vertices_all_key) {
        if (set_cardinality < cardinality[j][new_v]) {  // change parent(s) to key
          for (QueryVertexID existing_v : to_intersect_vertices[j]) {
            if (!(nxt_cover_node.cover_bits >> existing_v & 1)) {
              nxt_cover_node.cover_bits |= 1ULL << existing_v;
              nxt_cover_node.cover.push_back(existing_v);
            }
          }
        } else {  // change the new_v to key
          nxt_cover_node.cover_bits |= 1ULL << new_v;
          nxt_cover_node.cover.push_back(new_v);
        }
      }

      bool existing = false;
      for (const auto& cover_node : covers_[j]) {
        if (cover_node.cover_bits == nxt_cover_node.cover_bits) {
          existing = true;
          break;
        }
      }
      if (!existing) {
        covers_[j].emplace_back(nxt_cover_node);
      }
    }

    CoverNode last_cover_node = new_cover_node;
    for (uint32_t j = i; j > 0; --j) {
      QueryVertexID delete_v = matching_order_[j];
      if (!(last_cover_node.cover_bits >> delete_v & 1)) {
        for (QueryVertexID existing_v : to_intersect_vertices[j]) {
          const auto nbrs = query_graph_->getOutNeighbors(existing_v);
          uint32_t all_key = 1;
          for (uint32_t nbr_i = 0; nbr_i < nbrs.second; ++nbr_i) {
            QueryVertexID nbr_v = nbrs.first[nbr_i];
            if (existing_vertices[j].find(nbr_v) != existing_vertices[j].end()) {
              all_key &= (last_cover_node.cover_bits >> nbr_v) & 1;
            }
          }
          if (all_key) {
            last_cover_node.cover_bits -= 1ULL << existing_v;
            for (uint32_t it = 0; it < last_cover_node.cover.size(); ++it) {
              if (last_cover_node.cover[it] == existing_v) {
                last_cover_node.cover.erase(last_cover_node.cover.begin() + it);
                break;
              }
            }
          }
        }
      } else {
        last_cover_node.cover_bits -= 1ULL << delete_v;
        for (uint32_t it = 0; it < last_cover_node.cover.size(); ++it) {
          if (last_cover_node.cover[it] == delete_v) {
            last_cover_node.cover.erase(last_cover_node.cover.begin() + it);
            break;
          }
        }
      }
      bool existing = false;
      for (const auto& cover_node : covers_[j - 1]) {
        if (cover_node.cover_bits == last_cover_node.cover_bits) {
          existing = true;
          break;
        }
      }
      if (!existing) {
        covers_[j - 1].emplace_back(last_cover_node);
      }
    }
  }

  CoverNode first_level_node;
  bool key_existing = false;
  bool set_existing = false;
  for (const auto& cover_node : covers_[0]) {
    if (cover_node.cover_bits == 0) {
      set_existing = true;
    } else {
      CHECK_EQ(cover_node.cover_bits, 1ULL << matching_order_[0]);
      key_existing = true;
    }
  }
  // first node is set
  if (!set_existing) {
    first_level_node.cover_bits = 0;
    covers_[0].emplace_back(first_level_node);
  }
  // first node is key
  if (!key_existing) {
    first_level_node.cover_bits = 1ULL << matching_order_[0];
    first_level_node.cover.emplace_back(matching_order_[0]);
    covers_[0].emplace_back(first_level_node);
  }

  // populate cover nodes' parents
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    uint32_t v = 0;
    for (auto& cover_node : covers_[i]) {
      for (uint32_t j = 0; j < covers_[i - 1].size(); ++j) {
        const auto& last_level_cover_node = covers_[i - 1][j];
        if ((cover_node.cover_bits | last_level_cover_node.cover_bits) == cover_node.cover_bits) {
          cover_node.parents.emplace_back(j);
        }
      }
    }
  }
  // check whether there is a path from level 0 to level n
  std::vector<std::vector<uint64_t>> can_go(covers_.size());
  can_go[0].emplace_back(1);
  can_go[0].emplace_back(1);
  for (uint32_t i = 1; i < covers_.size(); ++i) {
    can_go[i].resize(covers_[i].size());
    std::fill(can_go[i].begin(), can_go[i].end(), 0);
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      auto& cover_node = covers_[i][j];
      for (uint32_t par : cover_node.parents) {
        can_go[i][j] += can_go[i - 1][par];
      }
    }
  }
  const uint32_t last = covers_.size() - 1;
  uint64_t sum = 0;
  for (uint32_t i = 0; i < covers_[last].size(); ++i) {
    sum += can_go[last][i];
  }
  CHECK_GT(sum, 0);
}

std::tuple<uint32_t, uint32_t, uint32_t> NaivePlanner::analyzeDynamicCoreCoverMWVC(
    const std::vector<QueryVertexID>& use_order) {
  // if any of the candidate cardinality is zero, there is no matching
  if (!hasValidCandidate()) {
    return std::tuple(0, 0, 0);
  }
  if (!use_order.empty()) {
    matching_order_ = use_order;
  }
  CHECK(!matching_order_.empty());

  // meters
  uint32_t sum_delayed_steps = 0;
  uint32_t sum_key_sizes = 0;
  uint32_t sum_key_to_set_changes = 0;

  auto subquery_vertices = matching_order_;
  auto candidate_cardinality = *candidate_cardinality_;

  // the static cover for the original query
  auto log_cardinality = logCardinality();
  WeightedBnB vertex_cover_solver(query_graph_, log_cardinality);
  vertex_cover_solver.computeVertexCover();
  auto& covers = vertex_cover_solver.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto choice = BnB::getSmallestCover(covers);
  std::vector<bool> current_cover_assignment;
  current_cover_assignment.resize(matching_order_.size());
  std::stringstream ss;
  for (QueryVertexID i = 0; i < current_cover_assignment.size(); ++i) {
    current_cover_assignment[i] = (covers[choice.first.front()][i] == 1);
    ss << ' ' << (covers[choice.first.front()][matching_order_[i]] == 1);
  }
  LOG(INFO) << "cover" << ss.str() << " from " << choice.first.size();
  sum_key_sizes += choice.second;
  // compute the key sizes for static cover
  if (current_cover_assignment[matching_order_.front()]) {
    sum_delayed_steps += matching_order_.size() - 1;  // the first vertex appears from step 1
  }
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    auto v = matching_order_[i];
    if (current_cover_assignment[v]) {
      sum_delayed_steps += (matching_order_.size() - i);  // v is key for all the rest steps
    }
  }

  // use mwvc to decide the dynamic cover for each subquery i
  for (uint32_t i = matching_order_.size() - 2; i > 0; --i) {
    subquery_vertices.resize(i + 1);
    candidate_cardinality.resize(i + 1);
    // compute a vertex cover of subquery i
    auto subquery = query_graph_->getInducedSubgraph(subquery_vertices);
    WeightedBnB vc_solver(&subquery, candidate_cardinality);
    vc_solver.computeVertexCover();
    auto& covers = vc_solver.getBestCovers();
    CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
    auto choice = BnB::getSmallestCover(covers);
    auto select_cover = covers[choice.first.front()];
    sum_key_sizes += choice.second;
    // for existing vertices 0 ~ i, check the change in key/set assignment
    std::stringstream ss;
    for (QueryVertexID j = 0; j <= i; ++j) {
      if (select_cover[j] == 1) {  // now assign j to key
        // add to sum_key_to_set_changes if change from key to set
        sum_key_to_set_changes += (!current_cover_assignment[matching_order_[j]]);
      }
      ss << ' ' << (select_cover[j] == 1);
      current_cover_assignment[matching_order_[j]] = (select_cover[j] == 1);
    }
    LOG(INFO) << "cover" << ss.str() << " from " << choice.first.size();
  }

  sum_delayed_steps -= sum_key_sizes;
  return std::tuple(sum_key_sizes, sum_delayed_steps, sum_key_to_set_changes);
}

}  // namespace circinus

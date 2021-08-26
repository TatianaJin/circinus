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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/k_core.h"
#include "algorithms/minimum_weight_vertex_cover.h"
#include "algorithms/vertex_cover.h"
#include "graph/query_graph.h"
#include "ops/order_generator.h"
#include "plan/execution_plan.h"
#include "utils/hashmap.h"

namespace circinus {

bool NaivePlanner::hasValidCandidate() const {
  uint32_t i = 0;
  for (auto& cardinality : candidate_cardinality_) {
    if (cardinality < 0) {
      LOG(INFO) << "candidate for vertex is empty " << i;
      return false;
    }
    ++i;
  }
  return true;
}

ExecutionPlan* NaivePlanner::generatePlan() {
  if (matching_order_.empty()) {
    return nullptr;
  }

  // get a smallest minimum weight vertex cover of the query graph
  auto log_cardinality = logCardinality();
  WeightedBnB vertex_cover_solver(query_graph_, log_cardinality);
  vertex_cover_solver.computeVertexCover();
  auto& covers = vertex_cover_solver.getBestCovers();
  CHECK_GT(covers.size(), 0);  // at least one cover should be obtained
  auto select_cover = covers[BnB::getSmallestCover(covers).first.front()];

  uint64_t cover_bits = 0;
  for (QueryVertexID i = 0; i < select_cover.size(); ++i) {
    cover_bits |= ((uint32_t)select_cover[i]) << i;
  }

  if (shortPlannerLog()) {
    std::stringstream ss;
    for (QueryVertexID i = 0; i < select_cover.size(); ++i) {
      ss << ' ' << select_cover[i];
    }
    LOG(INFO) << "Static cover table" << ss.str();
  }

  plan_.setQueryCoverBits(cover_bits);
  plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover);
  return &plan_;
}

ExecutionPlan* NaivePlanner::generatePlanWithoutCompression() {
  if (matching_order_.empty()) {
    return nullptr;
  }

  // all query vertices are in cover, so that no compression is done
  auto n_qvs = query_graph_->getNumVertices();
  std::vector<int> select_cover(n_qvs, 1);
  uint64_t cover_bits = 0;
  for (QueryVertexID i = 0; i < n_qvs; ++i) {
    cover_bits |= 1ULL << i;
  }

  plan_.setQueryCoverBits(cover_bits);
  plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover);
  return &plan_;
}

unordered_map<VertexID, double> NaivePlanner::dfsComputeCost(QueryVertexID qid, const uint64_t cover_bits,
                                                             std::vector<bool>& visited, std::vector<bool>& visited_cc,
                                                             const GraphBase* data_graph,
                                                             const std::vector<CandidateSetView>& candidate_views,
                                                             const std::vector<QueryVertexID>& cc,
                                                             bool with_traversal) {
  visited[qid] = true;
  visited_cc[cc[qid]] = true;
  unordered_map<VertexID, double> ret;
  for (VertexID v : candidate_views[qid]) {
    ret[v] = 1;
  }

  auto[qnbrs, qcnt] = query_graph_->getOutNeighbors(qid);
  if (!with_traversal) {
    /* For a DFS tree from qid in the given cc:
     * 1) initialize the weights of leaf candidates as 1;
     * 2) for each candidate, take sum of the weights of the neighbor candidates in each child tree node, and take
     * product of the sums;
     * 3) in case that a neighbor is a non-cover query vertex, each sum is averaged over the number of the neighbor
     * candidates assuming the neighbor candidates' neighborhood overlaps.
     */
    for (uint32_t i = 0; i < qcnt; ++i) {  // for each unvisited neighbor in the cc, recursively compute cost
      QueryVertexID nbr = qnbrs[i];
      // dfs into each cc from other different cc only one time
      if (visited[nbr] || cc[nbr] == DUMMY_QUERY_VERTEX || (cc[nbr] != cc[qid] && visited_cc[cc[nbr]])) {
        continue;
      }

      auto child_ret =
          dfsComputeCost(nbr, cover_bits, visited, visited_cc, data_graph, candidate_views, cc, with_traversal);
      for (auto& pair : ret) {  // for each candidate of qid, compute the weight from neighbor candidates
        double sum = 0;
        double num = 0;
        auto[nbrs, cnt] = data_graph->getOutNeighbors(pair.first);
        for (uint32_t j = 0; j < cnt; ++j) {
          auto pos = child_ret.find(nbrs[j]);
          if (pos != child_ret.end()) {
            sum += pos->second;
            num += 1;
          }
        }
        if (!(cover_bits >> nbr & 1) && num != 0) {
          sum /= num / 2;
        }
        pair.second *= sum;
      }
    }
  } else {
    /* For a DFS tree from qid in the given cc:
     * 1) initialize the weights of leaf candidates as 1;
     * 2) for each candidate at a tree node, compute the subtree cardinality at each child tree node, and take product;
     *   a) for a child node in cover, take sum of the weights of the neighbors in the child candidates;
     *   b) for a child node not in cover, traverse two hops and find the two-hop neighbor candidates, for each
     *      descendant node (which must be in cover) take sum of the neighbor candidate weights, and take product;
     * 3) traversal to child nodes in the cover is prioritized at each node to make sure nodes in the original
     * cover-only ccs are connected in the dfs tree without passing through a non-cover node.
     */
    for (uint32_t i = 0; i < qcnt; ++i) {
      QueryVertexID nbr = qnbrs[i];
      if (visited[nbr] || cc[nbr] == DUMMY_QUERY_VERTEX) {
        continue;
      }
      if (cover_bits >> nbr & 1) {
        auto child_ret =
            dfsComputeCost(nbr, cover_bits, visited, visited_cc, data_graph, candidate_views, cc, with_traversal);
        for (auto& pair : ret) {
          double sum = 0;
          auto[nbrs, cnt] = data_graph->getOutNeighbors(pair.first);
          for (uint32_t j = 0; j < cnt; ++j) {
            auto pos = child_ret.find(nbrs[j]);
            if (pos != child_ret.end()) {
              sum += pos->second;
            }
          }
          pair.second *= sum;
        }
      }
    }
    for (uint32_t i = 0; i < qcnt; ++i) {
      QueryVertexID nbr = qnbrs[i];
      if (visited[nbr] || cc[nbr] == DUMMY_QUERY_VERTEX) {
        continue;
      }
      if ((cover_bits >> nbr & 1) == 0) {
        visited[nbr] = true;
        auto[set_nbrs, set_cnt] = query_graph_->getOutNeighbors(nbr);
        unordered_set<VertexID> set_candidate(candidate_views[nbr].begin(), candidate_views[nbr].end());
        for (uint32_t j = 0; j < set_cnt; ++j) {  // for each two-hop descendant node, compute cost for each candidate
          QueryVertexID set_nbr = set_nbrs[j];
          if (visited[set_nbr] || cc[set_nbr] == DUMMY_QUERY_VERTEX) {
            continue;
          }
          auto child_ret =
              dfsComputeCost(set_nbr, cover_bits, visited, visited_cc, data_graph, candidate_views, cc, with_traversal);

          // two-hop traversal to find the cover node candidates
          for (auto& pair : ret) {
            unordered_set<VertexID> has_added;
            double sum = 0;
            auto[nbrs, cnt] = data_graph->getOutNeighbors(pair.first);
            for (uint32_t k = 0; k < cnt; ++k) {
              if (set_candidate.count(nbrs[k]) == 1) {
                auto[second_nbrs, second_cnt] = data_graph->getOutNeighbors(nbrs[k]);
                for (uint32_t l = 0; l < second_cnt; ++l) {
                  VertexID second_level_vid = second_nbrs[l];
                  auto pos = child_ret.find(second_level_vid);
                  if (pos != child_ret.end() && has_added.insert(second_level_vid).second) {
                    sum += pos->second;
                  }
                }
              }
            }
            pair.second *= sum;
          }
        }
      }
    }
  }
  return ret;
}

void NaivePlanner::getCoverCC(QueryVertexID qid, QueryVertexID cc_id, const uint64_t cover_bits,
                              std::vector<QueryVertexID>& cc) const {
  cc[qid] = cc_id;
  auto[qnbrs, qcnt] = query_graph_->getOutNeighbors(qid);
  for (uint32_t i = 0; i < qcnt; ++i) {
    QueryVertexID nbr = qnbrs[i];
    if ((cover_bits >> nbr & 1) == 0 || cc[nbr] != DUMMY_QUERY_VERTEX) {
      continue;
    }
    getCoverCC(nbr, cc_id, cover_bits, cc);
  }
}

void NaivePlanner::getMinimalConectedSubgraphWithAllKeys(std::vector<QueryVertexID>& cc, uint32_t level,
                                                         uint64_t cover_bits) {
  std::stringstream ss;
  if (verbosePlannerLog()) {  // debug log
    for (uint32_t i = 0; i <= level; ++i) {
      ss << matching_order_[i] << ", ";
    }
    DLOG(INFO) << "subquery vertices: " << ss.str();
  }
  uint32_t cc_count = 0;
  /* find connected components in cover */
  for (uint32_t i = 0; i <= level; ++i) {
    auto qid = matching_order_[i];
    if ((cover_bits >> qid & 1) == 0 || cc[qid] != DUMMY_QUERY_VERTEX) {
      continue;
    }
    ++cc_count;
    getCoverCC(qid, qid, cover_bits, cc);
  }
  if (verbosePlannerLog()) {  // debug log
    ss.str("");
    unordered_set<uint32_t> connected_components;
    for (uint32_t i = 0; i < matching_order_.size(); ++i) {
      if (cc[i] != DUMMY_QUERY_VERTEX) {
        connected_components.insert(cc[i]);
        ss << i << ':' << cc[i] << ", ";
      }
    }
    CHECK_EQ(connected_components.size(), cc_count);
    DLOG(INFO) << "keys: " << ss.str();
  }

  if (cc_count == 1) return;

  /* if there are multiple CCs, greedily pick non-cover vertices with smallest candidate cardinality to connect CCs */
  std::vector<uint32_t> increasing_idx(level + 1, 0);
  std::iota(increasing_idx.begin(), increasing_idx.end(), 0);
  sort(increasing_idx.begin(), increasing_idx.end(), [&](uint32_t l, uint32_t r) {
    return candidate_cardinality_[matching_order_[l]] < candidate_cardinality_[matching_order_[r]];
  });
  unordered_set<QueryVertexID> connecting_cover_cc;
  for (uint32_t i : increasing_idx) {
    QueryVertexID qid = matching_order_[i];
    if ((cover_bits >> qid & 1) == 0) {
      auto[nbr, cnt] = query_graph_->getOutNeighbors(qid);
      uint32_t current_cc_size = connecting_cover_cc.size();
      for (uint32_t j = 0; j < cnt; ++j) {
        if (cover_bits >> nbr[j] & 1) {
          CHECK_NE(cc[nbr[j]], DUMMY_QUERY_VERTEX) << nbr[j];
          connecting_cover_cc.insert(cc[nbr[j]]);
        }
      }
      if (connecting_cover_cc.size() > current_cc_size) {
        cc[qid] = qid;
        if (cc_count == connecting_cover_cc.size()) break;
      }
    }
  }

  if (verbosePlannerLog()) {  // debug log
    ss.str("");
    for (uint32_t i = 0; i < matching_order_.size(); ++i) {
      if (cc[i] != DUMMY_QUERY_VERTEX) {
        ss << i << ':' << cc[i] << ", ";
      }
    }
    DLOG(INFO) << "after adding set: " << ss.str();
  }
}

double NaivePlanner::estimateCardinality(const GraphBase* data_graph,
                                         const std::vector<CandidateSetView>& candidate_views, uint64_t cover_bits,
                                         uint32_t level) {
  double ret = 1;
  std::vector<bool> visited(matching_order_.size(), false);
  std::vector<bool> visited_cc(matching_order_.size(), false);
  std::vector<QueryVertexID> cc(matching_order_.size(), DUMMY_QUERY_VERTEX);
  getMinimalConectedSubgraphWithAllKeys(cc, level, cover_bits);

  for (uint32_t i = level + 1; i < matching_order_.size(); ++i) {
    visited[matching_order_[i]] = true;
  }

  /* dfs from a cover vertex to compute estimation */
  // TODO(tatiana): find one cover vertex with least candidates?
  for (auto qid : matching_order_) {
    if ((cover_bits >> qid & 1) == 0 || visited[qid]) {
      continue;
    }
    ret = 0;
    auto current_candidate_cars =
        dfsComputeCost(qid, cover_bits, visited, visited_cc, data_graph, candidate_views, cc, use_two_hop_traversal_);
    for (auto& candidate_car : current_candidate_cars) {
      ret += candidate_car.second;
    }
    break;
  }

  // DLOG(INFO) << "cost: " << ret;
  return ret;
}

std::vector<double> NaivePlanner::getParallelizingQueryVertexWeights(
    QueryVertexID partition_qv, const GraphBase* data_graph, const std::vector<CandidateSetView>& candidate_views,
    uint64_t cover_bits) {
  std::vector<bool> visited(matching_order_.size(), false);
  std::vector<bool> visited_cc(matching_order_.size(), false);
  std::vector<QueryVertexID> cc(matching_order_.size(), DUMMY_QUERY_VERTEX);
  // LOG(INFO) << cover_bits;
  getMinimalConectedSubgraphWithAllKeys(cc, matching_order_.size() - 1, cover_bits);

  // use partition qv as the starting vertex.
  //  1) partition qv is in the connected subgraph;
  //  2) partition qv is not in the connected subgraph:
  //     there must be a neighbor of pqv in the connected subgraph
  auto current_candidate_cars = dfsComputeCost(partition_qv, cover_bits, visited, visited_cc, data_graph,
                                               candidate_views, cc, use_two_hop_traversal_);

  std::vector<double> weights(candidate_views[partition_qv].size(), 0.0);
  uint32_t idx = 0;
  double sum = 0;
  for (VertexID v : candidate_views[partition_qv]) {
    auto pos = current_candidate_cars.find(v);
    DCHECK(pos != current_candidate_cars.end());
    weights[idx] = pos->second;
    sum += weights[idx];
    ++idx;
  }
  DLOG(INFO) << "parallelizing vertex weight sum " << sum;
  return std::move(weights);
}

double NaivePlanner::estimateExpandCost(const GraphBase* data_graph,
                                        const std::vector<CandidateSetView>* candidate_views,
                                        const std::vector<std::vector<double>>& car,
                                        const unordered_set<QueryVertexID>& existing_vertices, uint32_t level,
                                        uint32_t idx, uint32_t parent) {
  if (candidate_views == nullptr) {  // update the estimation of computation for set to key / expand into
    return car[level - 1][parent];
  }

  auto cover_bits = covers_[level][idx].cover_bits;
  auto parent_cover_bits = covers_[level - 1][parent].cover_bits;
  auto target_vertex = matching_order_[level];
  bool target_in_cover = cover_bits >> target_vertex & 1;
  uint32_t key_parent_cnt = 0, set_parent_cnt = 0;

  auto[nbrs, cnt] = query_graph_->getOutNeighbors(target_vertex);
  QueryVertexID set_parent = DUMMY_QUERY_VERTEX;
  for (uint32_t k = 0; k < cnt; ++k) {
    if (existing_vertices.count(nbrs[k]) != 0) {
      if ((parent_cover_bits >> nbrs[k] & 1) == 0) {
        if (set_parent == DUMMY_QUERY_VERTEX ||
            (*candidate_views)[nbrs[k]].size() < (*candidate_views)[set_parent].size()) {
          set_parent = nbrs[k];
        }
        ++set_parent_cnt;
      } else {
        ++key_parent_cnt;
      }
    }
  }

  if (FLAGS_intersection_count_coefficient) {  // compute the cost by intersection count
    // expand from key
    // if (key_parent_cnt == 0) key_parent_cnt = 1;
    double cost = key_parent_cnt * car[level - 1][parent];
    if (target_in_cover) {
      if (set_parent_cnt == 0) {  // key to key only
        return cost;
      }
      // computation cost for expand-into / expand-from-set
      if (key_parent_cnt == 0) {  // set to key
        auto key_bits = parent_cover_bits | (1 << set_parent);
        cost = estimateCardinality(data_graph, *candidate_views, key_bits, level - 1);
      }
      auto key_bits = parent_cover_bits | (1 << target_vertex);  // including target
      return cost + set_parent_cnt * estimateCardinality(data_graph, *candidate_views, key_bits, level);
    }
    if (set_parent_cnt > 0) {  // enumerate key expand to set
      // computation cost to expand from enumerated parents
      return cost + set_parent_cnt * estimateCardinality(data_graph, *candidate_views, cover_bits, level - 1);
    }
    return cost;  // key to set
  }

  if (target_in_cover) {                                       // compute the cost by the number of compressed groups
    auto key_bits = parent_cover_bits | (1 << target_vertex);  // including target
    return estimateCardinality(data_graph, *candidate_views, key_bits, level);
  } else if (set_parent_cnt > 0) {
    return estimateCardinality(data_graph, *candidate_views, cover_bits, level - 1);
  }
  return car[level - 1][parent];
}

ExecutionPlan* NaivePlanner::generatePlanWithDynamicCover(const GraphBase* data_graph,
                                                          const std::vector<CandidateSetView>* candidate_views) {
  if (verbosePlannerLog()) {
    std::stringstream ss;
    for (auto& view : *candidate_views) {
      ss << " " << view.size();
    }
    LOG(INFO) << "generatePlanWithDynamicCover candidates" << ss.str();
    LOG(INFO) << "use_two_hop_traversal_ " << use_two_hop_traversal_;
  }

  /* generate pruned dynamic cover search space */
  std::vector<std::vector<double>> cardinality_per_level(candidate_cardinality_.size(), candidate_cardinality_);
  generateCoverNode(cardinality_per_level);

  if (verbosePlannerLog()) {
    logCoverSpace();
  }

  std::vector<std::vector<double>> costs_car(covers_.size());  // dp table: cost of each subquery given a cover
  std::vector<std::vector<uint32_t>> pre(covers_.size());      // best parent of cover for backtracing
  std::vector<std::vector<double>> car(covers_.size());  // cardinality estimation of compressed groups for each step

  /* estimate the compressed group cardinality for each subquery given each cover */
  for (uint32_t i = 0; i < covers_.size(); ++i) {
    car[i].resize(covers_[i].size());
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      car[i][j] = estimateCardinality(candidate_cardinality_, data_graph, candidate_views, covers_[i][j], i);
    }
  }

  /* dynamic programming to compute costs for sequences */
  unordered_set<QueryVertexID> existing_vertices;
  existing_vertices.insert(matching_order_[0]);
  costs_car[0].resize(covers_[0].size(), 0);  // initialization
  for (uint32_t i = 1; i < covers_.size(); ++i) {
    costs_car[i].resize(covers_[i].size(), -1);
    pre[i].resize(covers_[i].size());
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      auto& cover_node = covers_[i][j];
      for (uint32_t par : cover_node.parents) {
        if (costs_car[i - 1][par] >= 0) {
          double cost = estimateExpandCost(data_graph, candidate_views, car, existing_vertices, i, j, par);
          if (costs_car[i][j] < 0 || costs_car[i - 1][par] + cost < costs_car[i][j]) {
            costs_car[i][j] = costs_car[i - 1][par] + cost;
            pre[i][j] = par;
          }
        }
      }
    }
    existing_vertices.insert(matching_order_[i]);
  }

  uint32_t last = covers_.size() - 1;
  /* select the best sequence and trace back */
  double mini_cost = -1;
  int best_idx = -1;
  for (uint32_t i = 0; i < covers_[last].size(); ++i) {
    if (costs_car[last][i] >= 0 && (mini_cost < 0 || mini_cost > costs_car[last][i])) {
      mini_cost = costs_car[last][i];
      best_idx = i;
    }
  }
  CHECK_NE(best_idx, -1) << "ERROR: can not get best plan index.";
  if (verbosePlannerLog()) {  // debug log
    LOG(INFO) << "===== Covers for query graph =====";
    for (uint32_t i = 0; i < covers_[last].size(); ++i) {
      LOG(INFO) << i << " cover " << covers_[last][i].getCoverTableString(matching_order_.size()) << " costs_car "
                << costs_car[last][i];
    }
    LOG(INFO) << "best last level cover idx " << best_idx << " cost " << costs_car[last][best_idx];
  }
  // best_idx = 5 < covers_[last].size() - 1 ? 5 : covers_.size() - 1;

  auto& select_cover_node = covers_[last][best_idx];
  auto select_cover = select_cover_node.getCoverTable(matching_order_.size());

  auto[level_become_key, step_costs] = traceBackCoverPath(costs_car, pre, best_idx);

  plan_.setStepCosts(std::move(step_costs));
  plan_.setQueryCoverBits(select_cover_node.cover_bits);
  plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover, level_become_key);
  return &plan_;
}

QueryVertexID NaivePlanner::selectParallelizingQueryVertex(
    uint64_t cover_bits, const std::vector<QueryVertexID>& parallelizing_qv_candidates) {
  for (auto v : parallelizing_qv_candidates) {  // pick a parallelizing qv from the candidates which is in cover
    if (cover_bits >> v & 1) {
      return v;
    }
  }
  // if none of the candidate is in cover, pick the first cover vertex in the order
  for (auto v : matching_order_) {
    if (cover_bits >> v & 1) {
      return v;
    }
  }
  CHECK(false) << "cannot find a valid parallelizing qv";
  return DUMMY_QUERY_VERTEX;
}

ExecutionPlan* NaivePlanner::generatePlanWithSampleExecution(const std::vector<std::vector<double>>& cardinality,
                                                             const std::vector<double>& level_cost) {
  CHECK_GT(covers_.size(), 0);
  std::vector<std::vector<double>> costs_car_sample_cost(covers_.size());
  std::vector<std::vector<double>> costs_car(covers_.size());
  std::vector<std::vector<uint32_t>> pre(covers_.size());

  /* estimate the compressed group cardinality for each subquery given each cover */
  std::vector<std::vector<double>> car(covers_.size());
  for (uint32_t i = 0; i < covers_.size(); ++i) {
    car[i].resize(covers_[i].size());
    for (uint32_t j = 0; j < covers_[i].size(); ++j) {
      car[i][j] = estimateCardinality(cardinality[i], nullptr, nullptr, covers_[i][j], i);
    }
  }

  /* dynamic programming to compute costs for sequences */
  costs_car_sample_cost[0].resize(covers_[0].size(), 0);  // initialization
  costs_car[0].resize(covers_[0].size(), 0);              // initialization
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

  const uint32_t last = covers_.size() - 1;
  /* select the best sequence and trace back */
  double mini_cost = -1;
  int best_idx = -1;
  for (uint32_t i = 0; i < covers_[last].size(); ++i) {
    LOG(INFO) << i << " cover " << covers_[last][i].getCoverTableString(last + 1) << " cardinality " << car[last][i]
              << ", costs_car_sample_cost " << costs_car_sample_cost[last][i] << ", costs_car " << costs_car[last][i];

    if (costs_car_sample_cost[last][i] >= 0 && (mini_cost < 0 || mini_cost > costs_car_sample_cost[last][i])) {
      mini_cost = costs_car_sample_cost[last][i];
      best_idx = i;
    }
  }
  CHECK_NE(best_idx, -1) << "ERROR: can not get best plan index.";
  LOG(INFO) << "best last level cover idx " << best_idx << "  " << car[last][best_idx] << " min cost " << mini_cost;
  auto select_cover = covers_[last][best_idx].getCoverTable(matching_order_.size());

  // trace back from the last level
  auto[level_become_key, step_costs] = traceBackCoverPath(costs_car, pre, best_idx);

  plan_.setStepCosts(std::move(step_costs));
  plan_.setQueryCoverBits(covers_[last][best_idx].cover_bits);
  plan_.populatePhysicalPlan(query_graph_, matching_order_, select_cover, level_become_key);

  return &plan_;
}

std::pair<unordered_map<QueryVertexID, uint32_t>, std::vector<double>> NaivePlanner::traceBackCoverPath(
    const std::vector<std::vector<double>>& costs_car, const std::vector<std::vector<uint32_t>>& pre, int best_idx) {
  unordered_map<QueryVertexID, uint32_t> level_become_key;
  std::vector<double> step_costs(matching_order_.size(), 0);

  // trace back from the last level
  uint32_t last = covers_.size() - 1;
  std::vector<uint32_t> best_path;
  best_path.reserve(covers_.size());
  best_path.emplace_back(best_idx);

  if (verbosePlannerLog()) {
    for (uint32_t idx = 0; idx < covers_[last].size(); ++idx) {  // debug log
      LOG(INFO) << "---------- sequence " << idx << " ----------";
      uint32_t last_idx = idx;
      for (uint32_t i = 1; i < matching_order_.size(); ++i) {
        LOG(INFO) << "[ " << covers_[last][last_idx].getCoverTableString(matching_order_.size()) << " ] "
                  << costs_car[last][last_idx] << " target " << matching_order_[i];
        last_idx = pre[last--][last_idx];
      }
      last = matching_order_.size() - 1;
    }
  }

  if (verbosePlannerLog()) {  // debug log
    LOG(INFO) << "---------- best sequence ----------";
  }
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    if (verbosePlannerLog()) {  // debug log
      LOG(INFO) << "[ " << covers_[last][best_idx].getCoverTableString(matching_order_.size()) << " ] "
                << costs_car[last][best_idx] << " target " << matching_order_[i];
    }
    best_idx = pre[last--][best_idx];
    best_path.emplace_back(best_idx);
  }

  std::reverse(best_path.begin(), best_path.end());
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    for (uint32_t j = 0; j <= i; ++j) {
      auto vid = matching_order_[j];
      if (covers_[i][best_path[i]].cover_bits >> vid & 1) {
        if (vid != matching_order_[i] && !(covers_[i - 1][best_path[i - 1]].cover_bits >> vid & 1)) {
          // if the first vertex becomes key in the next subquery, make it as key initially
          level_become_key.insert({vid, (i == 1) ? 0 : i});
          DLOG(INFO) << "add " << vid << " to key at " << i;
        }
        if (vid == matching_order_[i]) {
          level_become_key.insert({vid, i});
          DLOG(INFO) << "add " << vid << " to key at " << i;
        }
      }
    }
  }
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {  // the cost of each step along best_path
    step_costs[i] = costs_car[i][best_path[i]] - costs_car[i - 1][best_path[i - 1]];
    if (step_costs[i] == 0) {
      LOG(WARNING) << "step " << i << " - " << matching_order_[i]
                   << " cost is estimated as 0 = " << costs_car[i][best_path[i]] << " - "
                   << costs_car[i - 1][best_path[i - 1]] << ", "
                   << covers_[i][best_path[i]].getCoverTableString(matching_order_.size()) << " from "
                   << covers_[i - 1][best_path[i - 1]].getCoverTableString(matching_order_.size())
                   << " level become key "
                   << (level_become_key.count(matching_order_[i]) ? (int)level_become_key.at(matching_order_[i]) : -1);
    }
  }
  return {level_become_key, step_costs};
}

const std::vector<QueryVertexID>& NaivePlanner::generateOrder(const GraphBase* data_graph,
                                                              const GraphMetadata& metadata,
                                                              const std::vector<CandidateSetView>* candidate_views,
                                                              const std::vector<VertexID>& candidate_cardinality,
                                                              OrderStrategy order_strategy,
                                                              const std::vector<QueryVertexID>* use_order) {
  matching_order_.clear();
  if (hasValidCandidate()) {
    if (use_order == nullptr) {
      auto order_generator =
          OrderGenerator(data_graph, metadata, query_graph_, *candidate_views, candidate_cardinality);
      matching_order_ = order_generator.getOrder(order_strategy);
    } else {
      CHECK_EQ(use_order->size(), query_graph_->getNumVertices());
      matching_order_ = *use_order;
    }
  }
  return matching_order_;
}

std::vector<std::vector<int>> NaivePlanner::generateAnchorCovers(const std::vector<QueryVertexID>& subquery_vertices,
                                                                 const std::vector<double>& cardinality) {
  auto subquery = query_graph_->getInducedSubgraph(subquery_vertices);
  WeightedBnB vc_solver(&subquery, cardinality);
  vc_solver.computeVertexCover();
  return std::move(vc_solver.getBestCovers());
}

std::vector<std::vector<int>> NaivePlanner::generateAnchorCoversSettingParentAsKey(
    const std::vector<QueryVertexID>& subquery_vertices, const std::vector<double>& cardinality,
    const std::vector<uint32_t>& order_index, const std::vector<std::vector<QueryVertexID>>& parent_vertices_per_level,
    const std::vector<unordered_set<QueryVertexID>>& existing_vertices) {
  auto subquery_size = subquery_vertices.size();
  // no anchor cover for the full query due to no further expansion
  if (subquery_size == matching_order_.size()) {
    return {};
  }
  auto subquery = query_graph_->getInducedSubgraph(subquery_vertices);
  WeightedBnB vc_solver(&subquery, cardinality);
  // make parent query vertices at level subquery_size as keys, virtually removing them from subquery
  std::vector<QueryVertexID> pre_set_to_keys;
  // pre set to intersect vertices to key
  for (QueryVertexID to_intersect_vertex : parent_vertices_per_level[subquery_size]) {
    pre_set_to_keys.push_back(order_index[to_intersect_vertex]);
  }
  vc_solver.computeVertexCover(&pre_set_to_keys);
  auto& covers = vc_solver.getBestCovers();
  auto choice = BnB::getSmallestCover(covers);
  auto select_cover = covers[choice.first.front()];

  // set to_intersect_key to set if they can be set
  for (QueryVertexID to_intersect_vertex : parent_vertices_per_level[subquery_size]) {
    auto nbrs = query_graph_->getOutNeighbors(to_intersect_vertex);
    bool can_be_set = 1;
    for (uint32_t j = 0; j < nbrs.second; ++j) {
      QueryVertexID nbr = nbrs.first[j];
      if (existing_vertices[subquery_size].find(nbr) != existing_vertices[subquery_size].end()) {
        can_be_set &= select_cover[order_index[nbr]] == 1;
      }
    }
    if (can_be_set) {
      select_cover[order_index[to_intersect_vertex]] = 0;
    }
  }
  return {std::move(select_cover)};
}

void NaivePlanner::generateCoverNode(const std::vector<std::vector<double>>& cardinality) {
  if (!hasValidCandidate()) {
    return;
  }

  auto subquery_vertices = matching_order_;
  covers_.resize(matching_order_.size());

  // compute existing vertices for each level j when target is matching_order_[j]
  std::vector<unordered_set<QueryVertexID>> existing_vertices;
  existing_vertices.resize(matching_order_.size());
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    for (uint32_t j = i + 1; j < matching_order_.size(); ++j) {
      existing_vertices[j].insert(matching_order_[i]);
    }
  }

  // compute parent vertices for each level i when target is matching_order_[i]
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
  }

  // take log of the cardinality for computing weighted vertex cover
  // TODO(tatiana): more efficient implementation by reusing the prefix values? or the cardinality varies among
  // subqueries?
  std::vector<std::vector<double>> candidate_cardinality_by_match_order(matching_order_.size());
  std::vector<uint32_t> order_index(matching_order_.size(), 0);  // query vertex to matching order index
  for (uint32_t i = 0; i < matching_order_.size(); ++i) {
    order_index[matching_order_[i]] = i;
    candidate_cardinality_by_match_order[i].resize(i + 1);
    for (uint32_t j = 0; j < i + 1; ++j) {
      candidate_cardinality_by_match_order[i][j] = log2(cardinality[i][matching_order_[j]]);
    }
  }

  /* populate each level in the search space covers_ */
  for (uint32_t i = matching_order_.size() - 1; i > 0; --i) {
    subquery_vertices.resize(i + 1);
    // compute anchor vertex covers of subquery i
    auto covers = generateAnchorCovers(subquery_vertices, candidate_cardinality_by_match_order[i]);
    // auto covers = generateAnchorCoversSettingParentAsKey(subquery_vertices, candidate_cardinality_by_match_order[i],
    // order_index, to_intersect_vertices, existing_vertices);
    CHECK_GT(covers.size(), 0) << "level " << i << " doesn't have any cover.";  // at least one cover should be obtained
    for (const auto& select_cover : covers) {
      /* add cover to search space if not existing */
      CoverNode new_cover_node;
      for (uint32_t j = 0; j < select_cover.size(); ++j) {
        if (select_cover[j] == 1) {
          QueryVertexID v = matching_order_[j];
          new_cover_node.cover_bits |= 1ULL << v;
        }
      }

      if (verbosePlannerLog() && i == matching_order_.size() - 1) {  // debug log
        std::string s = "";
        std::vector<uint32_t> mapped_cover(matching_order_.size(), 0);
        for (uint32_t j = 0; j < select_cover.size(); ++j) {
          if (select_cover[j] == 1) {
            CHECK_LT(matching_order_[j], mapped_cover.size()) << matching_order_[j] << " " << mapped_cover.size();
            mapped_cover[matching_order_[j]] = 1;
          }
        }
        for (uint32_t j = 0; j < mapped_cover.size(); ++j) {
          s += std::to_string(mapped_cover[j]) + ",";
        }
        DLOG(INFO) << "subquery " << i << " cover " << s;
      }

      DLOG(INFO) << "------------------- level = " << i;

      addCover(new_cover_node, i);

      /* populate compatible covers for larger subqueries for new added cover node */
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
              nxt_cover_node.cover_bits |= 1ULL << existing_v;
            }
          } else {  // change the new_v to key
            nxt_cover_node.cover_bits |= 1ULL << new_v;
          }
        }

        addCover(nxt_cover_node, j);
      }

      /* populate compatible covers for smaller subqueries from new added cover node */
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
              CHECK_GE(last_cover_node.cover_bits, (1ULL << existing_v)) << i << " " << existing_v;
              last_cover_node.cover_bits -= 1ULL << existing_v;
            }
          }
        } else {
          CHECK_GE(last_cover_node.cover_bits, (1ULL << delete_v));
          last_cover_node.cover_bits -= 1ULL << delete_v;
        }
        addCover(last_cover_node, j - 1);
      }
    }
  }

  /* populate covers for the single-vertex subquery */
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
    covers_[0].emplace_back(first_level_node);
  }

  /* populate cover nodes' parents */
  for (uint32_t i = 1; i < matching_order_.size(); ++i) {
    for (auto& cover_node : covers_[i]) {
      for (uint32_t j = 0; j < covers_[i - 1].size(); ++j) {
        const auto& last_level_cover_node = covers_[i - 1][j];
        if ((cover_node.cover_bits | last_level_cover_node.cover_bits) == cover_node.cover_bits) {  // key containment
          // avoid enumerate key expand to key
          if ((cover_node.cover_bits >> matching_order_[i] & 1) &&
              (last_level_cover_node.cover_bits | (1 << matching_order_[i])) != cover_node.cover_bits) {
            continue;
          }
          if ((cover_node.cover_bits >> matching_order_[i] & 1) == 0) {
            auto mask = last_level_cover_node.cover_bits;
            for (auto x : to_intersect_vertices[i]) {
              mask |= 1 << x;
            }
            if (mask != cover_node.cover_bits) {
              continue;
            }
          }
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

void NaivePlanner::logCoverSpace() {
  for (uint32_t i = 0; i < covers_.size(); ++i) {
    std::stringstream ss;
    for (auto& cover : covers_[i]) {
      ss << '\t' << cover.getCoverTableString(matching_order_.size()) << " [";
      for (auto parent : cover.parents) {
        ss << ' ' << parent;
      }
      ss << "]";
    }
    LOG(INFO) << ss.str();
  }
}
}  // namespace circinus

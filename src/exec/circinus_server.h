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

#include <chrono>
#include <deque>
#include <fstream>
#include <string>

#include "glog/logging.h"

#include "exec/executor_manager.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "plan/planner.h"
#include "utils/hashmap.h"
#include "utils/query_utils.h"

namespace circinus {

struct QueryState {
  QueryContext query_context;
  Planner* planner = nullptr;

  QueryState(QueryGraph&& q, QueryConfig&& config, Graph* g) : query_context(std::move(q), std::move(config), g) {}
};

// TODO(tatiana): use a message queue?
class CircinusServer {
 protected:
  unordered_map<std::string, Graph> data_graphs_;
  std::deque<QueryState> active_queries_;
  std::vector<uint32_t> reusable_indices_;

  ExecutorManager executor_manager_;

 public:
  inline bool hasActiveQuery() const { return active_queries_.size() > reusable_indices_.size(); }
  inline uint32_t numDataGraphs() const { return data_graphs_.size(); }

  /** Load graph from binary file.
   * @param graph_path The binary file path .
   * @param name The name by which the loaded graph is referred to.
   */
  double loadGraphFromBinary(const std::string& graph_path, const std::string& name) {
    std::ifstream input(graph_path, std::ios::binary);
    CHECK(input.is_open());
    auto start_loading = std::chrono::steady_clock::now();
    Graph& data_graph = data_graphs_[name];
    data_graph.loadCompressed(input);
    auto end_loading = std::chrono::steady_clock::now();
    return ((double)std::chrono::duration_cast<std::chrono::microseconds>(end_loading - start_loading).count()) / 1e9;
  }

  /** Handles new Query, used together with prepareQuery(uint32_t).
   * @param graph_name The graph queried.
   * @param query_file The file storing the query graph.
   * @param query_config_str The string of the query config.
   * @returns The index of the new query.
   */
  uint32_t newQuery(const std::string& graph_name, const std::string& query_file, const std::string& query_config_str) {
    if (reusable_indices_.empty()) {
      active_queries_.emplace_back(QueryGraph(query_file), QueryConfig(query_config_str), &data_graphs_.at(graph_name));
      return active_queries_.size() - 1;
    }
    // reuse deleted index
    auto idx = reusable_indices_.back();
    active_queries_[idx].query_context =
        QueryContext(QueryGraph(query_file), QueryConfig(query_config_str), &data_graphs_.at(graph_name));
    reusable_indices_.pop_back();
    return idx;
  }

  /** Invokes Phase 1 planning and execution of a query.
   * @param query_index The index of the query to run.
   */
  void prepareQuery(uint32_t query_index) {
    auto& query_state = active_queries_[query_index];
    query_state.planner = new Planner(query_state.query_context);
    // phase 1: preprocessing
    auto plan = query_state.planner->generateCandidatePruningPlan();
    // asynchronous execution, the callback executeQuery() will be invoked when preprocessing finish
    executor_manager_.computeCandidates(query_index, plan);
  }

  /** Invokes Phase 2 planning and execution of a query.
   * @param query_index The index of the query to run.
   * @param candidate_sets The candidate sets of query vertices.
   */
  void executeQuery(uint32_t query_index, std::vector<std::vector<VertexID>>&& candidate_sets) {
    auto& query_state = active_queries_[query_index];
    DCHECK(query_state.planner != nullptr);
    // phase 2: query execution plan
    auto plan = query_state.planner->generateExecutionPlan();
    // asynchronous execution, the callback finishQuery() will be invoked when execution finish
    executor_manager_.executeQuery(query_index, plan);
  }

  /** Handles query results.
   */
  void finishQuery(uint32_t query_index, void* result) {
    delete active_queries_[query_index].planner;
    active_queries_[query_index].planner = nullptr;
    reusable_indices_.push_back(query_index);
    // TODO(tatiana): call back to client
  }
};

}  // namespace circinus

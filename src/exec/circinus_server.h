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
#include <memory>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "glog/logging.h"
#include "zmq.hpp"

#include "exec/candidate_pruning_plan_driver.h"
#include "exec/executor_manager.h"
#include "exec/threadsafe_queue.h"
#include "graph/graph.h"
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "plan/planner.h"
#include "utils/file_utils.h"
#include "utils/hashmap.h"
#include "utils/query_utils.h"

namespace circinus {

struct QueryState {
  QueryContext query_context;
  Planner* planner = nullptr;
  std::string client_addr;

  QueryState(QueryGraph&& q, QueryConfig&& config, Graph* g, GraphMetadata* m)
      : query_context(std::move(q), std::move(config), g, m) {}
};

class CircinusServer {
 protected:
  using Event = ServerEvent;

  unordered_map<std::string, std::pair<Graph, GraphMetadata>> data_graphs_;
  std::deque<QueryState> active_queries_;
  std::vector<uint32_t> reusable_indices_;

  ExecutorManager executor_manager_;

  ThreadsafeQueue<Event> event_queue_;
  std::thread client_server_thread_;

  zmq::context_t zmq_ctx_;
  // TODO(tatiana): clear disconnected sockets
  std::unordered_map<std::string, zmq::socket_t> sockets_to_clients_;

 public:
  CircinusServer() : executor_manager_(&event_queue_) {}
  ~CircinusServer() {
    LOG(INFO) << "CircinusServer shuts down now.";
    if (client_server_thread_.joinable()) {
      client_server_thread_.join();
      LOG(INFO) << "CircinusServer client port closed.";
    }
  }

  void Listen();
  void Serve(bool listen = true);
  auto& getZMQContext() { return zmq_ctx_; }

  inline void shutDown() { event_queue_.push(Event(Event::ShutDown)); }
  bool loadGraph(std::string&& gpath, std::string&& gname, std::string&& config = "", std::string&& client_addr = "");
  bool query(std::string&& game, std::string&& qpath, std::string&& config, std::string&& client_addr = "");

 protected:
  inline bool hasActiveQuery() const { return active_queries_.size() > reusable_indices_.size(); }
  inline uint32_t numDataGraphs() const { return data_graphs_.size(); }

  /* start of event handlers */

  /**
   * See newQuery(const std::string&, const std::string&, const std::string&) for values in event.args.
   */
  void handleNewQuery(const Event& event);

  /** Load data graph.
   * TODO(tatiana): support asynchronous parallel graph loading
   * @param event The event contains three args: graph path, graph name, and loading config.
   */
  void handleLoadGraph(const Event& event);

  void handleShutDown(const Event& event) {
    // TODO(tatiana): finish up before shutdown
  }

  void handleCandidatePhase(const Event& event);

  inline void handleExecutionPhase(const Event& event) {
    DCHECK(event.data != nullptr);
    finishQuery(event.query_id, event.data, "");
  }

  /* end of event handlers */

  /** Load graph from binary file.
   * @param graph_path The binary file path .
   * @param name The name by which the loaded graph is referred to.
   * @returns Graph loading time in seconds.
   */
  double loadGraphFromBinary(const std::string& graph_path, const std::string& name) {
    std::ifstream input(graph_path, std::ios::binary);
    if (!input.is_open()) {                // cannot open file
      if (Path::isRelative(graph_path)) {  // if relative path, try search under FLAGS_data_dir
        return loadGraphFromBinary(Path::join(FLAGS_data_dir, graph_path), name);
      }
      std::stringstream ss;
      ss << std::strerror(errno) << ": " << graph_path;
      throw std::runtime_error(ss.str());
    }
    auto start_loading = std::chrono::steady_clock::now();
    Graph data_graph;
    data_graph.loadUndirectedGraphFromBinary(input);
    data_graph.buildLabelIndex();  // TODO(tatiana): construct label index only when suitable
    auto memory = data_graph.getMemoryUsage();
    LOG(INFO) << "graph " << name << " takes " << memory.first / 1024 / 1024 << "MB";
    GraphMetadata meta(data_graph);  // TODO(tatiana): more statistics, support partitioned graph
    auto end_loading = std::chrono::steady_clock::now();
    data_graphs_.insert({name, std::make_pair(data_graph, meta)});
    return ((double)std::chrono::duration_cast<std::chrono::microseconds>(end_loading - start_loading).count()) / 1e6;
  }

  /** Handles new Query, used together with prepareQuery(uint32_t).
   * @param graph_name The graph queried.
   * @param query_file The file storing the query graph.
   * @param query_config_str The string of the query config.
   * @returns The index of the new query.
   */
  inline uint32_t newQuery(const std::string& graph_name, const std::string& query_file,
                           const std::string& query_config_str) {
    QueryGraph q(query_file);
    QueryConfig config(query_config_str);

    auto& graph = data_graphs_.at(graph_name);
    if (reusable_indices_.empty()) {
      active_queries_.emplace_back(std::move(q), std::move(config), &graph.first, &graph.second);
      return active_queries_.size() - 1;
    }
    // reuse deleted index
    auto idx = reusable_indices_.back();
    active_queries_[idx].query_context = QueryContext(std::move(q), std::move(config), &graph.first, &graph.second);
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
    // asynchronous execution, a CandidatePhase event will be generated when preprocessing finish
    executor_manager_.run(query_index, &query_state.query_context, std::make_unique<CandidatePruningPlanDriver>(plan));
  }

  /** Handles query results.
   */
  void finishQuery(uint32_t query_index, void* result, const std::string& error) {
    auto& query = active_queries_[query_index];
    delete query.planner;
    query.planner = nullptr;
    reusable_indices_.push_back(query_index);
    if (error.empty()) {
      if (query.query_context.query_config.output == "count") {
        replyToClient(query.client_addr, result, 8);
      } else {
        // TODO(tatiana): support other output option
        LOG(WARNING) << "Output option not implemented yet: " << query.query_context.query_config.output;
      }
    } else {
      replyToClient(query.client_addr, error);
    }
  }

  /**
   * Copy is incurred to copy data to zmq message buffer
   */
  bool replyToClient(const std::string& client_addr, const void* data, size_t size, bool success = true);
  inline bool replyToClient(const std::string& client_addr, const std::string& error_str) {
    if (!FLAGS_standalone || client_addr.empty()) LOG(WARNING) << error_str;
    return replyToClient(client_addr, error_str.data(), error_str.size(), false);
  }

  /**
   * Zero-copy reply where the passed data will be deleted by zmq socket
   */
  bool replyToClient(const std::string& client_addr, std::vector<char>* data);
};

}  // namespace circinus

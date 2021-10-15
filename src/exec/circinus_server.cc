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

#include "exec/circinus_server.h"

#include <memory>
#include <string>
#include <vector>

#include "zmq.hpp"
#include "zmq_addon.hpp"

#include "exec/execution_plan_driver.h"
#include "graph/partitioned_graph.h"
#include "utils/network_utils.h"

namespace circinus {

#define DEFAULT_PORT 55954
#define NUM_TRY_PORT 10

void CircinusServer::Listen() {
  client_server_thread_ = std::thread([this]() {
    zmq::socket_t socket(zmq_ctx_, ZMQ_PULL);
    // default port 55954, repeatedly try next if not available
    auto port = GetAvailablePort(DEFAULT_PORT, NUM_TRY_PORT);
    LOG(INFO) << "Circinus Server listening on port " << port;
    socket.bind("tcp://*:" + std::to_string(port));
    while (true) {
      zmq::multipart_t msg;
      msg.recv(socket);
      auto event_str = msg.popstr();
      if (event_str == "exit") {
        CHECK(msg.empty()) << msg.size();
        shutDown();
        break;
      } else if (event_str == "query") {
        DCHECK_GE(msg.size(), 3);
        auto graph_name = msg.popstr();
        auto query_path = msg.popstr();
        auto query_config = msg.popstr();
        if (msg.size() == 0) {
          query(std::move(graph_name), std::move(query_path), std::move(query_config));
        } else {  // need to reply to client
          query(std::move(graph_name), std::move(query_path), std::move(query_config), msg.popstr());
        }
      } else if (event_str == "load") {
        DCHECK_GE(msg.size(), 2);
        auto graph_path = msg.popstr();
        auto graph_name = msg.popstr();
        if (msg.size() == 0) {
          loadGraph(std::move(graph_path), std::move(graph_name));
        } else if (msg.size() == 1) {
          auto load_config = msg.popstr();
          loadGraph(std::move(graph_path), std::move(graph_name), std::move(load_config));
        } else {  // need to reply to client
          auto load_config = msg.popstr();
          auto client_addr = msg.popstr();
          loadGraph(std::move(graph_path), std::move(graph_name), std::move(load_config), std::move(client_addr));
        }
      }
    }
  });
}

void CircinusServer::Serve(bool listen) {
  if (listen) {
    Listen();
  }
  while (true) {
    auto event = event_queue_.waitAndPop();
    if (event.type == Event::ShutDown) {
      handleShutDown(event);
      break;
    }
    try {
      switch (event.type) {
      case Event::NewQuery: {
        handleNewQuery(event);
        break;
      }
      case Event::LoadGraph: {
        handleLoadGraph(event);
        break;
      }
      case Event::CandidatePhase: {
        handleCandidatePhase(event);
        break;
      }
      case Event::ExecutionPhase: {
        handleExecutionPhase(event);
        break;
      }
      case Event::TimeOut: {
        handleQueryTimeOut(event);
        break;
      }
      default:
        LOG(FATAL) << "Unknown event type " << event;
      }
    } catch (const std::exception& e) {
      if (event.client_addr.empty()) {
        LOG(WARNING) << e.what();
      } else {
        replyToClient(event.client_addr, std::string(e.what()));
      }
    }
  }
}

bool CircinusServer::loadGraph(std::string&& graph_path, std::string&& graph_name, std::string&& load_config,
                               std::string&& client_addr) {
  Event event(Event::LoadGraph);
  LOG(INFO) << "load " << graph_path << " " << graph_name << " " << load_config;
  event.args.reserve(3);
  event.args.emplace_back(std::move(graph_path));
  event.args.emplace_back(std::move(graph_name));
  event.args.emplace_back(std::move(load_config));
  event.client_addr = std::move(client_addr);
  event_queue_.push(std::move(event));
  return true;
}

bool CircinusServer::query(std::string&& graph_name, std::string&& query_path, std::string&& query_config,
                           std::string&& client_addr) {
  LOG(INFO) << "query " << graph_name << " " << query_path << " " << query_config;
  Event event(Event::NewQuery);
  event.args.reserve(3);
  event.args.emplace_back(std::move(graph_name));
  event.args.emplace_back(std::move(query_path));
  event.args.emplace_back(std::move(query_config));
  event.client_addr = std::move(client_addr);
  event_queue_.push(std::move(event));
  return true;
}

void CircinusServer::handleNewQuery(const Event& event) {
  if (data_graphs_.count(event.args[0]) == 0) {
    std::stringstream error;
    error << "Graph '" << event.args[0] << "' does not exist";
    replyToClient(event.client_addr, error.str());
    return;
  }
  auto query_idx = newQuery(event.args[0], event.args[1], event.args[2]);
  active_queries_[query_idx].client_addr = event.client_addr;
  prepareQuery(query_idx);
}

void CircinusServer::handleLoadGraph(const Event& event) {
  double load_time = 0;
  if (event.args.back().empty()) {
    // Graph loading is blocking, the server does not handle new requests when loading the graph
    load_time = loadGraphFromBinary(event.args[0], event.args[1]);
  } else {
    LOG(WARNING) << "The graph loading option is not implemented yet";
  }
  if (!event.client_addr.empty()) {
    replyToClient(event.client_addr, &load_time, sizeof(load_time));
  }
}

void CircinusServer::handleCandidatePhase(const Event& event) {
  DCHECK(event.data != nullptr);
  auto planner = active_queries_[event.query_id].planner;
  DCHECK(planner != nullptr);
  auto result = (const CandidateResult*)event.data;
  auto plan = planner->updateCandidatePruningPlan(result);
  std::unique_ptr<PlanDriver> plan_driver = nullptr;
  if (plan->isFinished()) {
    auto now = std::chrono::high_resolution_clock::now();
    active_queries_[event.query_id].filter_time = toSeconds(active_queries_[event.query_id].filter_start_time, now);
    // backtracking phase
    LOG(INFO) << "Candidate generation finished. Start backtracking.";
    auto plan = planner->generateExecutionPlan(result, FLAGS_num_cores > 1);
    if (active_queries_[event.query_id].query_context.graph_metadata->numPartitions() == 1 ||
        !active_queries_[event.query_id].query_context.query_config.use_partitioned_graph) {
      plan_driver = std::make_unique<MatchingParallelExecutionPlanDriver>(plan);
    } else if (active_queries_[event.query_id].query_context.query_config.candidate_pruning_strategy ==
               CandidatePruningStrategy::Online) {
      plan_driver = std::make_unique<OnlineQueryExecutionPlanDriver>(plan);
    } else {
      plan_driver = std::make_unique<ExecutionPlanDriver>(plan, &zmq_ctx_);
    }
    active_queries_[event.query_id].plan_time = toSeconds(now, std::chrono::high_resolution_clock::now());
    LOG(INFO) << "Generated backtracking plan in " << active_queries_[event.query_id].plan_time << " seconds";
  }
  executor_manager_.run(event.query_id, &active_queries_[event.query_id].query_context, std::move(plan_driver));
}

bool CircinusServer::replyToClient(const std::string& client_addr, zmq::multipart_t& msg) {
  if (client_addr.empty()) return true;
  auto pos = sockets_to_clients_.find(client_addr);
  if (pos == sockets_to_clients_.end()) {
    pos = sockets_to_clients_.insert({client_addr, zmq::socket_t(zmq_ctx_, ZMQ_PUSH)}).first;
    pos->second.setsockopt(ZMQ_LINGER, 0);  // do not linger after socket is closed
    pos->second.connect(client_addr);
  }
  try {
    msg.send(pos->second);
    return true;
  } catch (zmq::error_t& e) {
    LOG(WARNING) << e.what();
  }
  return false;
}

bool CircinusServer::replyToClient(const std::string& client_addr, const void* data, size_t size, bool success) {
  zmq::multipart_t msg;
  msg.pushtyp<bool>(success);
  msg.addmem(data, size);
  return replyToClient(client_addr, msg);
}

bool CircinusServer::replyToClient(const std::string& client_addr, std::vector<char>* data) {
  zmq::multipart_t reply;
  reply.pushtyp<bool>(true);
  zmq::message_t msg(data->data(), data->size(), [](void*, void* hint) { delete (std::vector<char>*)hint; }, data);
  reply.add(std::move(msg));
  return replyToClient(client_addr, reply);
}

double CircinusServer::loadGraphFromBinary(const std::string& graph_path, const std::string& name) {
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
  auto data_graph = GraphBase::loadGraphFromBinary(input);
  auto graph = dynamic_cast<Graph*>(data_graph.get());
  GraphMetadata meta;
  if (graph != nullptr) {
    graph->buildLabelIndex();  // TODO(tatiana): construct label index only when suitable
    meta = GraphMetadata(*graph, true);
  } else {
    meta = GraphMetadata(*dynamic_cast<ReorderedPartitionedGraph*>(data_graph.get()));
  }
  auto end_loading = std::chrono::steady_clock::now();
  auto memory = data_graph->getMemoryUsage();
  LOG(INFO) << "graph " << name << " takes " << memory.first / 1024 / 1024 << "MB";
  data_graphs_.insert({name, std::make_pair(std::move(data_graph), std::move(meta))});
  return ((double)std::chrono::duration_cast<std::chrono::microseconds>(end_loading - start_loading).count()) / 1e6;
}

inline uint32_t CircinusServer::newQuery(const std::string& graph_name, const std::string& query_file,
                                         const std::string& query_config_str) {
  QueryGraph q(query_file);
  QueryConfig config(query_config_str);
  DLOG(INFO) << "Query execution mode " << ((uint16_t)config.mode);

  auto& graph = data_graphs_.at(graph_name);
  if (reusable_indices_.empty()) {
    active_queries_.emplace_back(std::move(q), std::move(config), graph.first.get(), &graph.second,
                                 std::chrono::steady_clock::now() + config.time_limit);
    return active_queries_.size() - 1;
  }
  // reuse deleted index
  auto idx = reusable_indices_.back();
  active_queries_[idx].query_context = QueryContext(std::move(q), std::move(config), graph.first.get(), &graph.second,
                                                    std::chrono::steady_clock::now() + config.time_limit);
  reusable_indices_.pop_back();
  return idx;
}

void CircinusServer::finishQuery(uint32_t query_index, void* result, const std::string& error) {
  auto& query = active_queries_[query_index];
  delete query.planner;
  query.planner = nullptr;
  reusable_indices_.push_back(query_index);
  if (error.empty()) {
    if (result == nullptr) {
      std::string reply = "no matching";
      replyToClient(query.client_addr, reply.data(), reply.size(), true);
    } else if (query.query_context.query_config.output == "count") {
      auto res = reinterpret_cast<QueryResult*>(result);
      res->filter_time = query.filter_time;
      res->plan_time = query.plan_time;
      auto reply = res->toString();
      replyToClient(query.client_addr, reply.data(), reply.size(), true);
    } else if (query.query_context.query_config.output == "plan") {
      auto& reply = *reinterpret_cast<std::string*>(result);
      replyToClient(query.client_addr, reply.data(), reply.size(), true);
    } else if (query.query_context.query_config.output == "profile_count") {
      auto& reply = *reinterpret_cast<QueryResult*>(result);
      reply.filter_time = query.filter_time;
      reply.plan_time = query.plan_time;
      zmq::multipart_t msg;
      msg.pushtyp<bool>(true);
      msg.addstr(reply.profile_info);
      msg.addstr(reply.toString());
      replyToClient(query.client_addr, msg);
    } else {
      // TODO(tatiana): support other output option
      LOG(WARNING) << "Output option not implemented yet: " << query.query_context.query_config.output;
    }
  } else {
    replyToClient(query.client_addr, error);
  }
  // make sure the pointer result is not needed before clearQuery
  executor_manager_.clearQuery(query_index, error);
}

void CircinusServer::consolidateConfigs(QueryContext& qctx) {
  auto& conf = qctx.query_config;

  // use online for seed-given labeled queries and unlabeled queries, use cfl otherwise
  if (conf.candidate_pruning_strategy == CandidatePruningStrategy::Adaptive) {
    if (conf.seed.first != DUMMY_QUERY_VERTEX || qctx.query_graph.isUnlabeled()) {
      conf.candidate_pruning_strategy = CandidatePruningStrategy::Online;
    } else {
      conf.candidate_pruning_strategy = CandidatePruningStrategy::CFL;
    }
  }

  if (conf.isProfileMode()) {
    conf.output = "profile_count";
  }

  // CandidatePruningStrategy::Online: only get candidates for start query vertex if seed is not given
  if (conf.candidate_pruning_strategy == CandidatePruningStrategy::Online) {
    // consolidate expand config
    FLAGS_candidate_set_intersection = 3;
    LOG(INFO) << "set FLAGS_candidate_set_intersection = 3 for online cps";
    // consolidate order strategy
    if (conf.order_strategy != OrderStrategy::Online && conf.order_strategy != OrderStrategy::None) {
      LOG(INFO) << "set order_strategy = OrderStrategy::Online for online cps";
      conf.order_strategy = OrderStrategy::Online;  // Degree sorting
    }
  }
}

void CircinusServer::prepareQuery(uint32_t query_index) {
  auto& query_state = active_queries_[query_index];
  query_state.planner = new Planner(query_state.query_context);
  consolidateConfigs(query_state.query_context);
  if (query_state.query_context.query_config.mode == QueryMode::Explain) {  // dry run
    auto cardinality = query_state.planner->estimateCardinality();
    // FIXME(by)
    // auto plan = query_state.planner->generateExecutionPlan(&cardinality);
    // auto str = plan->toString();
    // query_state.query_context.query_config.output = "plan";
    // finishQuery(query_index, &str, "");
  } else {  // enter actual execution phases
    // TODO(tatiana): refactor to merge both branches
    if (query_state.query_context.query_config.seed.first != DUMMY_QUERY_VERTEX &&
        query_state.query_context.query_config.order_strategy == OrderStrategy::Online) {
      LOG(INFO) << "Skipped candidate generation. Start backtracking.";

      auto now = std::chrono::high_resolution_clock::now();
      auto plan = query_state.planner->generateExecutionPlan(query_state.query_context.query_config.seed);
      query_state.plan_time = toSeconds(now, std::chrono::high_resolution_clock::now());
      LOG(INFO) << "Generated backtracking plan in " << query_state.plan_time << " seconds";

      if (plan == nullptr) {
        finishQuery(query_index, nullptr, "");
        return;
      }

      std::unique_ptr<PlanDriver> plan_driver = std::make_unique<OnlineQueryExecutionPlanDriver>(plan);
      executor_manager_.run(query_index, &query_state.query_context, std::move(plan_driver));
    } else {
      query_state.filter_start_time = std::chrono::high_resolution_clock::now();

      // phase 1: preprocessing
      CandidatePruningPlan* plan = query_state.planner->generateCandidatePruningPlan();
      if (plan->isFinished()) {  // no candidate generation TODO(tatiana): in case of seed data vertex
        auto plan = query_state.planner->generateExecutionPlan((CandidateResult*)nullptr, FLAGS_num_cores > 1);
        executor_manager_.run(query_index, &query_state.query_context,
                              std::make_unique<ExecutionPlanDriver>(plan, &zmq_ctx_));
        return;
      }
      // asynchronous execution, a CandidatePhase event will be generated when preprocessing finish
      executor_manager_.run(query_index, &query_state.query_context,
                            std::make_unique<CandidatePruningPlanDriver>(plan));
    }
  }
}

}  // namespace circinus

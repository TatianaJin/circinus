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
#include "utils/network_utils.h"

namespace circinus {

void CircinusServer::Listen() {
  client_server_thread_ = std::thread([this]() {
    zmq::socket_t socket(zmq_ctx_, ZMQ_PULL);
    // socket.bind("tcp://*:" + std::to_string(GetAvailablePort()));
    socket.bind("tcp://*:55954");
    while (true) {
      zmq::multipart_t msg;
      msg.recv(socket);
      auto event_str = msg.popstr();
      if (event_str == "exit") {
        CHECK(msg.empty());
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
        DCHECK_GE(msg.size(), 3);
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
  auto plan = planner->updateCandidatePruningPlan((const std::vector<VertexID>*)event.data);
  std::unique_ptr<PlanDriver> plan_driver = nullptr;
  if (plan->isFinished()) {
    // backtracking phase
    LOG(INFO) << "Candidate generation finished. Start backtracking.";
    auto plan = planner->generateExecutionPlan((const std::vector<VertexID>*)event.data);
    if (active_queries_[event.query_id].query_context.graph_metadata->numPartitions() == 1) {
      plan_driver = std::make_unique<MatchingParallelExecutionPlanDriver>(plan);
    } else {
      plan_driver = std::make_unique<ExecutionPlanDriver>(plan);
    }
  }
  executor_manager_.run(event.query_id, &active_queries_[event.query_id].query_context, std::move(plan_driver));
}

bool CircinusServer::replyToClient(const std::string& client_addr, const void* data, size_t size, bool success) {
  if (client_addr.empty()) return true;
  auto pos = sockets_to_clients_.find(client_addr);
  if (pos == sockets_to_clients_.end()) {
    pos = sockets_to_clients_.insert({client_addr, zmq::socket_t(zmq_ctx_, ZMQ_PUSH)}).first;
    pos->second.setsockopt(ZMQ_LINGER, 0);  // do not linger after socket is closed
    pos->second.connect(client_addr);
  }
  try {
    zmq::multipart_t msg;
    msg.pushtyp<bool>(success);
    msg.addmem(data, size);
    msg.send(pos->second);
    return true;
  } catch (zmq::error_t& e) {
    LOG(WARNING) << e.what();
  }
  return false;
}

bool CircinusServer::replyToClient(const std::string& client_addr, std::vector<char>* data) {
  if (client_addr.empty()) return true;
  auto pos = sockets_to_clients_.find(client_addr);
  if (pos == sockets_to_clients_.end()) {
    pos = sockets_to_clients_.insert({client_addr, zmq::socket_t(zmq_ctx_, ZMQ_PUSH)}).first;
    pos->second.setsockopt(ZMQ_LINGER, 0);  // do not linger after socket is closed
    pos->second.connect(client_addr);
  }
  try {
    zmq::multipart_t reply;
    reply.pushtyp<bool>(true);
    zmq::message_t msg(data->data(), data->size(), [](void*, void* hint) { delete (std::vector<char>*)hint; }, data);
    reply.add(std::move(msg));
    reply.send(pos->second);
    return true;
  } catch (zmq::error_t& e) {
    LOG(WARNING) << e.what();
  }
  return false;
}

}  // namespace circinus

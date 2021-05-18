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

#include <string_view>

#include "zmq.hpp"
#include "zmq_addon.hpp"

#include "utils/network_utils.h"

namespace circinus {

void CircinusServer::Listen() {
  client_server_thread_ = std::thread([this]() {
    zmq::socket_t socket(zmq_ctx_, ZMQ_PULL);
    socket.bind("tcp://*:" + std::to_string(GetAvailablePort()));
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
        if (msg.size() == 3) {
          query(msg.popstr(), msg.popstr(), msg.popstr());  // graph name, query path, query config
        } else {                                            // need to reply to client
          query(msg.popstr(), msg.popstr(), msg.popstr(), msg.popstr());
        }
      } else if (event_str == "load") {
        DCHECK_GE(msg.size(), 3);
        if (msg.size() == 3) {                                  // need to reply to client
          loadGraph(msg.popstr(), msg.popstr(), msg.popstr());  // graph path, graph name, load config
        } else {
          loadGraph(msg.popstr(), msg.popstr(), msg.popstr(), msg.popstr());
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
    auto event = event_queue_.WaitAndPop();
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
      default:
        LOG(FATAL) << "Unknown event type " << event;
      }
    } catch (const std::exception& e) {
      replyToClient(event.client_addr, std::string(e.what()));
    }
  }
}

bool CircinusServer::loadGraph(std::string&& graph_path, std::string&& graph_name, std::string&& load_config,
                               std::string&& client_addr) {
  Event event(Event::LoadGraph);
  event.args.reserve(3);
  event.args.emplace_back(std::move(graph_path));
  event.args.emplace_back(std::move(graph_name));
  event.args.emplace_back(std::move(load_config));
  event.client_addr = std::move(client_addr);
  event_queue_.Push(std::move(event));
  return true;
}

bool CircinusServer::query(std::string&& graph_name, std::string&& query_path, std::string&& query_config,
                           std::string&& client_addr) {
  Event event(Event::NewQuery);
  event.args.reserve(3);
  event.args.emplace_back(std::move(graph_name));
  event.args.emplace_back(std::move(query_path));
  event.args.emplace_back(std::move(query_config));
  event.client_addr = std::move(client_addr);
  event_queue_.Push(std::move(event));
  return true;
}

void CircinusServer::handleNewQuery(const Event& event) {
  if (data_graphs_.count(event.args[0]) == 0) {
    std::stringstream error;
    error << "Graph '" << event.args[0] << "'does not exist";
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

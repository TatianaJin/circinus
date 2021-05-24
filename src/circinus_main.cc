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

#include <string_view>

#include <iostream>

#include "glog/logging.h"
#include "zmq.hpp"
#include "zmq_addon.hpp"

#include "exec/circinus_server.h"
#include "utils/flags.h"

using circinus::CircinusServer;

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_log_dir == "") {
    google::LogToStderr();
  }

  CircinusServer server;
  if (FLAGS_standalone) {
    LOG(INFO) << "Circinus server starts in standalone mode";
    std::thread cli([&server]() {
      std::string tok;
      std::vector<std::string> args;
      zmq::socket_t sock(server.getZMQContext(), ZMQ_PULL);
      sock.bind("inproc://client");
      auto handler = [ args = &args, &server, &sock ]() {
        if (args->front() == "exit") {
          server.shutDown();
          return true;
        }
        if (args->front() == "load") {
          if (args->size() < 3) {
            LOG(WARNING) << "Wrong number of arguments for load: load <graph_path> <graph_name>";
          } else if (args->size() == 3) {
            server.loadGraph(std::move((*args)[1]), std::move((*args)[2]), "", "inproc://client");
          } else {
            server.loadGraph(std::move((*args)[1]), std::move((*args)[2]), std::move((*args)[3]), "inproc://client");
          }
          zmq::multipart_t reply;
          reply.recv(sock);
          auto success = reply.poptyp<bool>();
          if (success) {
            std::cout << "Loaded graph in " << reply.poptyp<double>() << " seconds" << std::endl;
          } else {
            auto msg = reply.pop();
            std::cerr << std::string_view((char*)msg.data(), msg.size()) << std::endl;
          }
        } else if (args->front() == "query") {
          if (args->size() != 4) {
            LOG(WARNING) << "Wrong number of arguments for query: query <graph_name> <query_path> <query_config>";
          } else {
            server.query(std::move((*args)[1]), std::move((*args)[2]), std::move((*args)[3]), "inproc://client");
          }
          zmq::multipart_t reply;
          reply.recv(sock);
          auto success = reply.poptyp<bool>();
          if (success) {
            // TODO(tatiana): time and result
            std::cout << "Query finished" << std::endl;
          } else {
            auto msg = reply.pop();
            std::cerr << std::string_view((char*)msg.data(), msg.size()) << std::endl;
          }
        } else {
          LOG(WARNING) << "Unknown command " << args->front();
          LOG(INFO) << "Available commands: exit | load | query";
        }
        args->clear();
        std::cerr << "Circinus> ";
        return false;
      };
      std::cerr << "Circinus> ";
      while (std::cin >> tok) {
        if (tok.front() == ';') {
          if (args.empty()) {  // empty statement
            std::cerr << "Circinus> ";
            continue;
          }
          if (handler()) break;
        } else if (tok.back() == ';') {
          args.emplace_back(tok.data(), tok.size() - 1);
          if (handler()) break;
        } else {
          args.push_back(std::move(tok));
        }
      }
    });
    server.Serve(false);
    cli.join();
  } else {
    server.Serve();
  }
  return 0;
}

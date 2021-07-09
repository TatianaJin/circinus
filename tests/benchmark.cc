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

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#ifdef WITH_GPERF
#include "gperftools/profiler.h"
#endif
#include "zmq.hpp"
#include "zmq_addon.hpp"

#include "exec/circinus_server.h"
#include "utils/file_utils.h"


using circinus::CircinusServer;
using circinus::QueryVertexID;


#define BATCH_SIZE FLAGS_batch_size
#define toSeconds(start, end) \
  (((double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1e9)

DEFINE_string(output_file, "", "The output file path");
/**
 * Datasets:
 * dblp, eu2005, hprd, human, patents, wordnet, yeast, youtube
 */
DEFINE_string(dataset, "dblp", "The dataset to use");
DEFINE_string(query_mode, "dense", "Dense or sparse query");
DEFINE_uint64(query_size, 8, "The query size.");
DEFINE_uint64(match_limit, 1e5, "The limit of matches to find");
DEFINE_uint64(query_index, 1, "The index of query in the same category");
DEFINE_string(match_order, "", "Matching order");
DEFINE_string(filter, "nlf", "filter");
DEFINE_string(profile_file, "", "profile file");
DEFINE_string(vertex_cover, "static", "Vertex cover strategy: static, dynamic, all");
DEFINE_string(batch_file, "", "Batch query file");
DEFINE_bool(bipartite_graph, false, "use bipartite graph or not");


class QueryConfig {
 public:
  std::string dataset;
  uint32_t query_size;
  std::string query_mode;
  uint32_t query_index;
  std::string match_order;

  bool skipConfig() const { return skip_; }

  void readNextConfig(std::istream& str) {
    skip_ = false;
    std::string::size_type pos = 0;
    auto old_pos = pos;
    std::getline(str, line_);
    // parse dataset
    if ((pos = line_.find(',', pos)) == std::string::npos) {
      skip_ = true;
      return;
    }
    dataset = std::string(&line_[old_pos], pos - old_pos);
    if (dataset.substr(0, 7).compare("dataset") == 0) {
      skip_ = true;
      return;
    }
    old_pos = ++pos;
    // parse query_size
    CHECK_NE((pos = line_.find(',', pos)), std::string::npos);
    line_[pos] = '\0';
    query_size = std::atoi(&line_[old_pos]);
    old_pos = ++pos;
    // parse query_mode
    CHECK_NE((pos = line_.find(',', pos)), std::string::npos);
    query_mode = std::string(&line_[old_pos], pos - old_pos);
    old_pos = ++pos;
    // parse query_index
    CHECK_NE((pos = line_.find(',', pos)), std::string::npos);
    line_[pos] = '\0';
    query_index = std::atoi(&line_[old_pos]);
    old_pos = ++pos;
    // parse match_order
    pos = line_.find(',', pos);
    if (pos == line_.npos) {
      pos = line_.size();
    }
    if (line_[pos - 1] == ' ') {
      match_order = std::string(&line_[old_pos], pos - old_pos - 1);
    } else {
      match_order = std::string(&line_[old_pos], pos - old_pos);
    }
  }

  std::ostream& toString(std::ostream& ss) {
    if (skip_) return ss;
    ss << "-dataset " << dataset << " -query_size " << query_size << " -query_mode " << query_mode << " -query_index "
       << query_index << " -match_order " << match_order << std::endl;
    return ss;
  }

  friend std::istream& operator>>(std::istream& str, QueryConfig& config) {
    config.readNextConfig(str);
    return str;
  }

 private:
  std::string line_;
  bool skip_ = false;
};

#define ADDRESS "inproc://benchmark"

class Benchmark {
  CircinusServer* server_;
  zmq::socket_t sock_;

 public:
  explicit Benchmark(CircinusServer& server) : server_(&server), sock_(server.getZMQContext(), ZMQ_PULL) {
    sock_.bind(ADDRESS);
  }

  static inline std::string getQueryPath(const std::string& dataset, uint32_t query_size, const std::string& query_mode,
                                         uint32_t index, std::ostream* out = nullptr) {
    if (out != nullptr) {
      (*out) << dataset << ',' << query_size << ',' << query_mode << ',' << index << ',';
    }
    return circinus::Path::join(
        FLAGS_data_dir, dataset, "query_graph",
        "query_" + query_mode + "_" + std::to_string(query_size) + "_" + std::to_string(index) + ".graph");
  }

  static inline std::string getGraphPath(const std::string& dataset) {
    return circinus::Path::join(FLAGS_data_dir, dataset, "data_graph", dataset + ".graph.bin");
  }

  void loadDataset(const std::string& dataset, double& load_time) {
    server_->loadGraph(getGraphPath(dataset), std::string(dataset), "", ADDRESS);
    zmq::multipart_t reply;
    reply.recv(sock_);
    auto success = reply.poptyp<bool>();
    if (success) {
      load_time = reply.poptyp<double>();
    } else {
      auto msg = reply.pop();
      LOG(FATAL) << std::string_view((char*)msg.data(), msg.size()) << std::endl;
    }
  }

  void run(const std::string& dataset, uint32_t query_size, const std::string& query_mode, uint32_t index,
           const std::string& match_order, std::ostream* out) {
    // get query path and log config
    auto query_path = getQueryPath(dataset, query_size, query_mode, index, out);
    // run query
    std::stringstream config;
    config << "cps=" << FLAGS_filter << ",cs=" << FLAGS_vertex_cover << ",limit=" << FLAGS_match_limit
           << ",mo=" << match_order;
    server_->query(std::string(dataset), std::move(query_path), config.str(), ADDRESS);
    // log result
    zmq::multipart_t reply;
    reply.recv(sock_);
    auto success = reply.poptyp<bool>();
    auto msg = reply.pop();
    if (success) {
      (*out) << std::string_view((char*)msg.data(), msg.size());
    } else {
      LOG(WARNING) << std::string_view((char*)msg.data(), msg.size());
    }
    (*out) << std::endl;
  }
};

void run_benchmark(const std::string& query_file, Benchmark& benchmark, std::ostream* out) {
  (*out) << "dataset,query_size,query_mode,query_index,elapsed_execution_time,filter_time,plan_time,"
            "enumerate_time,n_embeddings,"
            "order\n";
  std::ifstream query_f(FLAGS_batch_file);
  LOG(INFO) << "batch query file" << query_file;
  CHECK(query_f.is_open());
  QueryConfig config;
  std::string dataset = "";
  double load_time = 0;
  std::string profile_prefix = "/data/share/users/byli/circinus/evaluation/profile/";
  while (query_f >> config) {
    if (config.skipConfig()) continue;
    config.toString(LOG(INFO) << ">>>>>>>>>>>>>>>>> query config -vertex_cover " << FLAGS_vertex_cover << " -filter "
                              << FLAGS_filter << " -match_limit " << FLAGS_match_limit << ' ');
    // TODO(tatiana): change profile_file from system flag to a query config
    // load graph if not cached
    FLAGS_profile_file = profile_prefix + config.dataset + '_' + config.query_mode + '_' +
                         std::to_string(config.query_size) + '_' + std::to_string(config.query_index) + '_' +
                         FLAGS_filter;
    LOG(INFO) << "-------------" << FLAGS_profile_file;
    if (config.dataset != dataset) {
      dataset = config.dataset;
      benchmark.loadDataset(dataset, load_time);
    }
    benchmark.run(config.dataset, config.query_size, config.query_mode, config.query_index, config.match_order, out);
  }
}

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_log_dir.empty()) {
    google::LogToStderr();
  }

  // header: dataset,query_size,query_mode,query_index,load_time,filter_time,plan_time,enumerate_time,n_embeddings,order
  std::ostream* out;
  std::ofstream fstream;
  if (!FLAGS_output_file.empty()) {
    fstream.open(FLAGS_output_file, std::ios::app);
    CHECK(fstream.is_open());
    out = &fstream;
  } else {
    out = &std::cout;
  }
  if (FLAGS_match_limit == 0) {
    FLAGS_match_limit = ~0ull;
  }

  CircinusServer server;
  std::thread server_thread([&server]() { server.Serve(false); });

  Benchmark benchmark(server);
  if (!FLAGS_batch_file.empty()) {
    run_benchmark(FLAGS_batch_file, benchmark, out);
    server.shutDown();
    server_thread.join();
    return 0;
  }

  double load_time = 0;
  benchmark.loadDataset(FLAGS_dataset, load_time);
  benchmark.run(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, FLAGS_match_order, out);

  fstream.close();

  server.shutDown();
  server_thread.join();
  return 0;
}

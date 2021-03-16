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
#include <fstream>
#include <string>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "exec/thread_pool.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/expand_edge_operator.h"
#include "ops/filters.h"
#include "ops/scans.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::LDFScan;
using circinus::NaivePlanner;
using circinus::NLFFilter;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::Task;
using circinus::ThreadPool;
using circinus::VertexID;

#define BATCH_SIZE FLAGS_batch_size

DEFINE_string(data_dir, "/data/share/project/haxe/data/subgraph_matching_datasets", "The directory of datasets");
DEFINE_string(dataset, "dblp", "The dataset to use");
DEFINE_string(query_mode, "dense", "Dense or sparse query");
DEFINE_uint64(query_size, 8, "The query size.");
DEFINE_uint64(query_index, 1, "The index of query in the same category");
DEFINE_string(match_order, "", "Matching order");

class Benchmark {
 protected:
  const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                              "patents", "wordnet", "yeast", "youtube"};

 public:
  void run(const std::string& dataset, uint32_t query_size, const std::string& query_mode, uint32_t index) {
    auto graph_path = dataset + "/data_graph/" + dataset + ".graph";
    auto query_path = dataset + "/query_graph/query_" + query_mode + "_" + std::to_string(query_size) + "_" +
                      std::to_string(index) + ".graph";
    run(graph_path, query_path);
    std::cout << ',' << dataset << ',' << query_size << ',' << query_mode << ',' << index << std::endl;
  }

 protected:
  std::vector<std::vector<VertexID>> getCandidateSets(const Graph& g, const QueryGraph& q) {
    std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      candidates[v].reserve(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
      LDFScan scan(&q, v, &g);
      NLFFilter filter(&q, v);
      std::vector<VertexID> buffer;
      buffer.reserve(BATCH_SIZE);
      while (scan.Scan(&buffer, BATCH_SIZE) > 0) {
        filter.Filter(g, buffer, &candidates[v]);
        buffer.clear();
      }
    }
    return candidates;
  }

  std::vector<QueryVertexID> getOrder(const std::string& order_str, uint32_t size) {
    std::vector<QueryVertexID> order;
    if (order_str.empty()) return order;
    order.reserve(size);
    QueryVertexID v = 0;
    for (auto c : order_str) {
      if (c == ' ') {
        order.push_back(v);
        v = 0;
      } else {
        v = v * 10 + (c - '0');
      }
    }
    order.push_back(v);
    CHECK_EQ(order.size(), size);
    return order;
  }

  void run(const std::string& graph_path, const std::string& query_path) {
    Graph g(FLAGS_data_dir + "/" + graph_path);  // load data graph
    LOG(INFO) << "========================";
    LOG(INFO) << "graph " << graph_path << " query " << query_path;
    QueryGraph q(FLAGS_data_dir + "/" + query_path);  // load query graph
    auto use_order = getOrder(FLAGS_match_order, q.getNumVertices());
    auto candidates = getCandidateSets(g, q);  // get candidates for each query vertex
    std::vector<double> candidate_cardinality;
    candidate_cardinality.reserve(candidates.size());
    for (auto& set : candidates) {
      double size = set.size();
      // use log to transform product to sum for mwvc
      candidate_cardinality.push_back(std::log(size));
    }
    NaivePlanner planner(&q, &candidate_cardinality);
    auto delayed_steps = planner.analyzeDynamicCoreCoverEager(use_order);
    auto& order = planner.getMatchingOrder();
    auto dc = planner.analyzeDynamicCoreCoverMWVC();
    std::stringstream ss;
    for (auto v : order) {
      ss << v << ' ';
    }
    // order,delayed_steps_eager,key_sizes_mwvc,delayed_steps_mwvc,key_to_set_mwvc
    auto[key_sizes_mwvc, delayed_steps_mwvc, key_to_set_mwvc] = dc;
    std::cout << ss.str() << ',' << delayed_steps << ',' << key_sizes_mwvc << ',' << delayed_steps_mwvc << ','
              << key_to_set_mwvc;
  }
};

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_log_dir == "") {
    google::LogToStderr();
  }
  Benchmark benchmark;
  benchmark.run(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index);
  return 0;
}

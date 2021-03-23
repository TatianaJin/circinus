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
#include <string>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gperftools/profiler.h"
#include "gtest/gtest.h"

#include "exec/thread_pool.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/expand_edge_operator.h"
#include "ops/filters.h"
#include "ops/order.h"
#include "ops/scans.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/hashmap.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::LDFScan;
using circinus::NaivePlanner;
using circinus::NLFFilter;
using circinus::CFLFilter;
using circinus::CFLOrder;
using circinus::DPISOFilter;
using circinus::OrderBase;
using circinus::TSOOrder;
using circinus::TSOFilter;
using circinus::GQLFilter;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::Task;
using circinus::ThreadPool;
using circinus::VertexID;
using circinus::unordered_map;
using circinus::unordered_set;

#define BATCH_SIZE FLAGS_batch_size
#define toSeconds(start, end) \
  (((double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1e9)

DEFINE_string(data_dir, "/data/share/project/haxe/data/subgraph_matching_datasets", "The directory of datasets");
DEFINE_string(output_file, "", "The output file path");
DEFINE_string(dataset, "dblp", "The dataset to use");
DEFINE_string(query_mode, "dense", "Dense or sparse query");
DEFINE_uint64(query_size, 8, "The query size.");
DEFINE_uint64(match_limit, 1e5, "The limit of matches to find");
DEFINE_uint64(query_index, 1, "The index of query in the same category");
DEFINE_string(match_order, "", "Matching order");
DEFINE_string(filter, "nlf", "filter");

class Benchmark1 {
 protected:
  const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                              "patents", "wordnet", "yeast", "youtube"};

 public:
  void run(const std::string& dataset, uint32_t query_size, const std::string& query_mode, uint32_t index,
           std::ostream* out) {
    auto graph_path = dataset + "/data_graph/" + dataset + ".graph";
    auto query_path = dataset + "/query_graph/query_" + query_mode + "_" + std::to_string(query_size) + "_" +
                      std::to_string(index) + ".graph";

    (*out) << dataset << ',' << query_size << ',' << query_mode << ',' << index << ':';
    run(graph_path, query_path, out);
  }

 protected:
  std::vector<std::vector<VertexID>> getCandidateSets(const Graph& g, const QueryGraph& q) {
    std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
    std::vector<uint32_t> candidate_size(q.getNumVertices());
    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      candidates[v].reserve(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
      LDFScan scan(&q, v, &g);
      NLFFilter filter(&q, v);
      std::vector<VertexID> buffer;
      buffer.reserve(BATCH_SIZE);
      while (scan.Scan(&buffer, BATCH_SIZE) > 0) {
        if (FLAGS_filter == "dpiso") {
          candidates[v].insert(candidates[v].end(), buffer.begin(), buffer.end());
        } else
          filter.Filter(g, buffer, &candidates[v]);
        buffer.clear();
      }
      candidate_size[v] = candidates[v].size();
    }
    if (FLAGS_filter == "cfl") {
      CFLOrder cfl_order;
      QueryVertexID start_vertex = cfl_order.getStartVertex(&g, &q, candidate_size);
      LOG(INFO) << "cfl order get start vertex " << start_vertex;
      CFLFilter cfl_filter(&q, &g, start_vertex);
      cfl_filter.Filter(candidates);
    } else if (FLAGS_filter == "dpiso") {
      OrderBase dpiso_order;
      QueryVertexID start_vertex = dpiso_order.getStartVertex(&g, &q, candidate_size);
      LOG(INFO) << "dpiso order get start vertex " << start_vertex;
      DPISOFilter dpiso_filter(&q, &g, start_vertex);
      dpiso_filter.Filter(candidates);
    } else if (FLAGS_filter == "tso") {
      TSOOrder tso_order;
      QueryVertexID start_vertex = tso_order.getStartVertex(&g, &q, candidate_size);
      LOG(INFO) << "tso order get start vertex " << start_vertex;
      TSOFilter tso_filter(&q, &g, start_vertex);
      tso_filter.Filter(candidates);
    } else if (FLAGS_filter == "gql") {
      std::vector<std::vector<VertexID>> gql_candidates;
      unordered_map<QueryVertexID, unordered_set<VertexID>> gql_map;
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        unordered_set<VertexID> gql_set;
        for (auto i : candidates[v]) {
          gql_set.insert(i);
        }
        gql_map[v] = std::move(gql_set);
      }
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        GQLFilter gql_filter(&q, v, &gql_map);
        gql_filter.preFilter(g, candidates[v]);
      }
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        GQLFilter gql_filter(&q, v, &gql_map);
        std::vector<VertexID> tempvec;
        gql_filter.Filter(g, candidates[v], &tempvec);
        gql_candidates.push_back(std::move(tempvec));
      }
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        LOG(INFO) << "vertex " << v << " " << candidate_size[v] << "/" << gql_candidates[v].size();
      }
      return gql_candidates;
    }

    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      LOG(INFO) << "vertex " << v << " " << candidate_size[v] << "/" << candidates[v].size();
    }
    return candidates;
  }

  void run(const std::string& graph_path, const std::string& query_path, std::ostream* out) {
    auto start_loading = std::chrono::steady_clock::now();
    Graph g(FLAGS_data_dir + "/" + graph_path);  // load data graph
    auto end_loading = std::chrono::steady_clock::now();
    LOG(INFO) << "========================";
    LOG(INFO) << "graph " << graph_path << " query " << query_path;
    QueryGraph q(FLAGS_data_dir + "/" + query_path);  // load query graph
    auto candidates = getCandidateSets(g, q);         // get candidates for each query vertex
    std::stringstream ss;
    for (auto v : candidates) {
      ss << v.size() << ' ';
    }
    (*out) << ss.str() << '\n';
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

  Benchmark1 benchmark;
  // header: dataset,query_size,query_mode,query_index,load_time,filter_time,plan_time,enumerate_time,n_embeddings,order
  std::ostream* out;
  if (FLAGS_output_file != "") {
    std::ofstream fstream;
    fstream.open(FLAGS_output_file);
    CHECK(fstream.is_open());
    out = &fstream;
    benchmark.run(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
    fstream.close();
  } else {
    out = &std::cout;
    benchmark.run(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  }

  return 0;
}

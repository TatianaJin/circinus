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
#ifdef WITH_GPERF
#include "gperftools/profiler.h"
#endif
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

const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                            "patents", "wordnet", "yeast", "youtube"};
const char data_dir[] = "/data/share/project/haxe/data/subgraph_matching_datasets";
const char answer_dir[] = "/data/share/users/qlma/circinus-test/answer/";
const std::vector<int> query_size_list = {4, 8, 12, 16, 20, 24, 32};
const std::vector<std::string> query_mode_list = {"dense", "sparse"};
const std::pair<int, int> query_index_range = {1, 200};

std::vector<std::vector<VertexID>> getCandidateSets(const Graph& g, const QueryGraph& q,
                                                    const std::string& filter_str) {
  const char* filter_cstr = filter_str.c_str();
  std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
  std::vector<uint32_t> candidate_size(q.getNumVertices());
  for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
    candidates[v].reserve(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
    LDFScan scan(&q, v, &g);
    NLFFilter filter(&q, v);
    std::vector<VertexID> buffer;
    buffer.reserve(BATCH_SIZE);
    while (scan.Scan(&buffer, BATCH_SIZE) > 0) {
      if (strcmp(filter_cstr, "dpiso") == 0) {
        candidates[v].insert(candidates[v].end(), buffer.begin(), buffer.end());
      } else {
        filter.Filter(g, buffer, &candidates[v]);
      }
      buffer.clear();
    }
    candidate_size[v] = candidates[v].size();
  }
  if (strcmp(filter_cstr, "cfl") == 0) {
    CFLOrder cfl_order;
    QueryVertexID start_vertex = cfl_order.getStartVertex(&g, &q, candidate_size);
    // LOG(INFO) << "cfl order get start vertex " << start_vertex;
    CFLFilter cfl_filter(&q, &g, start_vertex);
    cfl_filter.Filter(candidates);
  } else if (strcmp(filter_cstr, "dpiso") == 0) {
    OrderBase dpiso_order;
    QueryVertexID start_vertex = dpiso_order.getStartVertex(&g, &q, candidate_size);
    // LOG(INFO) << "dpiso order get start vertex " << start_vertex;
    DPISOFilter dpiso_filter(&q, &g, start_vertex);
    dpiso_filter.Filter(candidates);
  } else if (strcmp(filter_cstr, "tso") == 0) {
    TSOOrder tso_order;
    QueryVertexID start_vertex = tso_order.getStartVertex(&g, &q, candidate_size);
    // LOG(INFO) << "tso order get start vertex " << start_vertex;
    TSOFilter tso_filter(&q, &g, start_vertex);
    tso_filter.Filter(candidates);
  } else if (strcmp(filter_cstr, "gql") == 0) {
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
    //   for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
    //     LOG(INFO) << "vertex " << v << " " << candidate_size[v] << "/" << gql_candidates[v].size();
    //}
    return gql_candidates;
  }

  //  for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
  //   LOG(INFO) << "vertex " << v << " " << candidate_size[v] << "/" << candidates[v].size();
  //}
  return candidates;
}
void run(const std::string& dataset, const std::string& filter, std::vector<std::string>& answers) {
  auto graph_path = dataset + "/data_graph/" + dataset + ".graph";
  auto data_dir_str = std::string(data_dir);
  Graph g(data_dir_str + "/" + graph_path);  // load data graph
                                         //  LOG(INFO) << "========================";
                                         // LOG(INFO) << "graph " << graph_path << " query " << query_path;
  int index = 0;
  for (auto& query_size : query_size_list)
    for (auto& query_mode : query_mode_list)
      for (int i = query_index_range.first; i <= query_index_range.second; ++i) {
        auto query_path = dataset + "/query_graph/query_" + query_mode + "_" + std::to_string(query_size) + "_" +
                          std::to_string(i) + ".graph";
        auto query_dir = data_dir_str + "/" + query_path;
        std::ifstream infile(query_dir);
        if (!infile) continue;
        QueryGraph q(query_dir);  // load query graph
        std::stringstream ss;
        ss << dataset << ',' << query_size << ',' << query_mode << ',' << i << ':';
        auto candidates = getCandidateSets(g, q, filter);  // get candidates for each query vertex
        for (auto v : candidates) {
          ss << v.size() << ' ';
        }
        std::string result_str = ss.str();
        result_str.pop_back();
        std::string expect_str = answers[index++];
        EXPECT_EQ(result_str, expect_str);
      }
}

bool getAnswers(const std::string& filter, const std::string& dataset, std::vector<std::string>& answers) {
  auto answer_path = std::string(answer_dir) + filter;
  std::ifstream in(answer_path);
  if (in) {
    std::string line;
    while (getline(in, line)) {
      if (line.rfind(dataset, 0) == 0) answers.emplace_back(line);
    }
    return true;
  } else {
    return false;
  }
}

void filterTest(std::string filter, std::string dataset) {
  std::vector<std::string> answers;
  bool res = getAnswers(filter, dataset, answers);
  EXPECT_EQ(res, true);
  run(dataset, filter, answers);
}

TEST(TestCFLFilterCandidates, dblp) { filterTest("cfl", "dblp"); }
TEST(TestCFLFilterCandidates, eu2005) { filterTest("cfl", "eu2005"); }
TEST(TestCFLFilterCandidates, hprd) { filterTest("cfl", "hprd"); }
TEST(TestCFLFilterCandidates, human) { filterTest("cfl", "human"); }
TEST(TestCFLFilterCandidates, patents) { filterTest("cfl", "patents"); }
TEST(TestCFLFilterCandidates, wordnet) { filterTest("cfl", "wordnet"); }
TEST(TestCFLFilterCandidates, yeast) { filterTest("cfl", "yeast"); }
TEST(TestCFLFilterCandidates, youtube) { filterTest("cfl", "youtube"); }
TEST(TestTSOFilterCandidates, dblp) { filterTest("tso", "dblp"); }
TEST(TestTSOFilterCandidates, eu2005) { filterTest("tso", "eu2005"); }
TEST(TestTSOFilterCandidates, hprd) { filterTest("tso", "hprd"); }
TEST(TestTSOFilterCandidates, human) { filterTest("tso", "human"); }
TEST(TestTSOFilterCandidates, patents) { filterTest("tso", "patents"); }
TEST(TestTSOFilterCandidates, wordnet) { filterTest("tso", "wordnet"); }
TEST(TestTSOFilterCandidates, yeast) { filterTest("tso", "yeast"); }
TEST(TestTSOFilterCandidates, youtube) { filterTest("tso", "youtube"); }
TEST(TestGQLFilterCandidates, dblp) { filterTest("gql", "dblp"); }
TEST(TestGQLFilterCandidates, eu2005) { filterTest("gql", "eu2005"); }
TEST(TestGQLFilterCandidates, hprd) { filterTest("gql", "hprd"); }
TEST(TestGQLFilterCandidates, human) { filterTest("gql", "human"); }
TEST(TestGQLFilterCandidates, patents) { filterTest("gql", "patents"); }
TEST(TestGQLFilterCandidates, wordnet) { filterTest("gql", "wordnet"); }
TEST(TestGQLFilterCandidates, yeast) { filterTest("gql", "yeast"); }
TEST(TestGQLFilterCandidates, youtube) { filterTest("gql", "youtube"); }
TEST(TestDPISOFilterCandidates, dblp) { filterTest("dpiso", "dblp"); }
TEST(TestDPISOFilterCandidates, eu2005) { filterTest("dpiso", "eu2005"); }
TEST(TestDPISOFilterCandidates, hprd) { filterTest("dpiso", "hprd"); }
TEST(TestDPISOFilterCandidates, human) { filterTest("dpiso", "human"); }
TEST(TestDPISOFilterCandidates, patents) { filterTest("dpiso", "patents"); }
TEST(TestDPISOFilterCandidates, wordnet) { filterTest("dpiso", "wordnet"); }
TEST(TestDPISOFilterCandidates, yeast) { filterTest("dpiso", "yeast"); }
TEST(TestDPISOFilterCandidates, youtube) { filterTest("dpiso", "youtube"); }

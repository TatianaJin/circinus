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
#include <random>
#include <sstream>
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
#include "ops/filters.h"
#include "ops/logical_filters.h"
#include "ops/operators.h"
#include "ops/order.h"
#include "ops/scans.h"
#include "ops/types.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/hashmap.h"
#include "utils/profiler.h"
#include "utils/utils.h"
#include "ops/filter_and_order.h"
using circinus::FilterAndOrder;

using circinus::CompressedSubgraphs;
using circinus::ExecutionConfig;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::GraphType;
using circinus::GraphMetadata;
using circinus::NaivePlanner;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::Task;
using circinus::ThreadPool;
using circinus::VertexID;
using circinus::Profiler;
using circinus::CoverNode;
using circinus::QueryType;
using circinus::TraverseOperator;
using circinus::INVALID_VERTEX_ID;
// logical filter
using circinus::LogicalCFLFilter;
using circinus::LogicalGQLFilter;
using circinus::LogicalNLFFilter;
using circinus::LogicalTSOFilter;
using circinus::LogicalDPISOFilter;
using circinus::LogicalNeighborhoodFilter;

// physical filter
using circinus::NeighborhoodFilter;
using circinus::NLFFilter;
using circinus::GQLFilter;

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
        FilterAndOrder fao(&g, &q,filter);
        auto candidates = fao.getCandidateSets();  // get candidates for each query vertex
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

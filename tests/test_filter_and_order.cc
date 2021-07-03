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
#include "ops/scans.h"
#include "ops/stateful_filter_and_order.h"
#include "ops/types.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/hashmap.h"
#include "utils/profiler.h"
#include "utils/utils.h"
using circinus::StatefulFilterAndOrder;

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
const char cs_answer_dir[] = "/data/share/users/qlma/circinus-test/answer/";
const char order_answer_dir[] = "/data/share/users/qlma/circinus-test/order/answer/";
const std::vector<int> query_size_list = {4, 8, 12, 16, 20, 24, 32};
const std::vector<std::string> query_mode_list = {"dense", "sparse"};
const std::pair<int, int> query_index_range = {1, 200};

void run(const std::string& dataset, const std::string& filter, std::vector<std::string>& CSAnswers,
         std::vector<std::string>& OAnswers) {
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
        std::stringstream ss1, ss2;
        ss1 << dataset << ',' << query_size << ',' << query_mode << ',' << i << ':';
        ss2 << dataset << ',' << query_size << ',' << query_mode << ',' << i << ':';
        StatefulFilterAndOrder fao(&g, &q, filter);
        auto candidates = fao.getCandidateSets();  // get candidates for each query vertex
        auto order = fao.getOrder();
        for (auto v : candidates) {
          ss1 << v.size() << ' ';
        }
        for (auto v : order) {
          ss2 << v << ' ';
        }
        std::string result_str1 = ss1.str(), result_str2 = ss2.str();
        result_str1.pop_back();
        result_str2.pop_back();
        std::string expect_str1 = CSAnswers[index];
        std::string expect_str2 = OAnswers[index];
        ++index;
        EXPECT_EQ(result_str1, expect_str1);
        EXPECT_EQ(result_str2, expect_str2);
      }
}
bool getAnswers(const std::string& path, const std::string& dataset, std::vector<std::string>& answers) {
  std::ifstream in(path);
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
bool getCSAnswers(const std::string& filter, const std::string& dataset, std::vector<std::string>& CSAnswers) {
  auto answer_path = std::string(cs_answer_dir) + filter;
  return getAnswers(answer_path, dataset, CSAnswers);
}
bool getOrderAnswers(const std::string& filter, const std::string& dataset, std::vector<std::string>& OAnswers) {
  auto answer_path = std::string(order_answer_dir) + filter;
  return getAnswers(answer_path, dataset, OAnswers);
}

void filterTest(std::string filter, std::string dataset) {
  std::vector<std::string> CSAnswers, OAnswers;
  bool res = getCSAnswers(filter, dataset, CSAnswers);
  EXPECT_EQ(res, true);
  res = getOrderAnswers(filter, dataset, OAnswers);
  EXPECT_EQ(res, true);
  run(dataset, filter, CSAnswers, OAnswers);
}

TEST(TestCFLFAO, dblp) { filterTest("cfl", "dblp"); }
TEST(TestCFLFAO, eu2005) { filterTest("cfl", "eu2005"); }
TEST(TestCFLFAO, hprd) { filterTest("cfl", "hprd"); }
TEST(TestCFLFAO, human) { filterTest("cfl", "human"); }
TEST(TestCFLFAO, patents) { filterTest("cfl", "patents"); }
TEST(TestCFLFAO, wordnet) { filterTest("cfl", "wordnet"); }
TEST(TestCFLFAO, yeast) { filterTest("cfl", "yeast"); }
TEST(TestCFLFAO, youtube) { filterTest("cfl", "youtube"); }
TEST(TestTSOFAO, dblp) { filterTest("tso", "dblp"); }
TEST(TestTSOFAO, eu2005) { filterTest("tso", "eu2005"); }
TEST(TestTSOFAO, hprd) { filterTest("tso", "hprd"); }
TEST(TestTSOFAO, human) { filterTest("tso", "human"); }
TEST(TestTSOFAO, patents) { filterTest("tso", "patents"); }
TEST(TestTSOFAO, wordnet) { filterTest("tso", "wordnet"); }
TEST(TestTSOFAO, yeast) { filterTest("tso", "yeast"); }
TEST(TestTSOFAO, youtube) { filterTest("tso", "youtube"); }
TEST(TestGQLFAO, dblp) { filterTest("gql", "dblp"); }
TEST(TestGQLFAO, eu2005) { filterTest("gql", "eu2005"); }
TEST(TestGQLFAO, hprd) { filterTest("gql", "hprd"); }
TEST(TestGQLFAO, human) { filterTest("gql", "human"); }
TEST(TestGQLFAO, patents) { filterTest("gql", "patents"); }
TEST(TestGQLFAO, wordnet) { filterTest("gql", "wordnet"); }
TEST(TestGQLFAO, yeast) { filterTest("gql", "yeast"); }
TEST(TestGQLFAO, youtube) { filterTest("gql", "youtube"); }
TEST(TestDPISOFAO, dblp) { filterTest("dpiso", "dblp"); }
TEST(TestDPISOFAO, eu2005) { filterTest("dpiso", "eu2005"); }
TEST(TestDPISOFAO, hprd) { filterTest("dpiso", "hprd"); }
TEST(TestDPISOFAO, human) { filterTest("dpiso", "human"); }
TEST(TestDPISOFAO, patents) { filterTest("dpiso", "patents"); }
TEST(TestDPISOFAO, wordnet) { filterTest("dpiso", "wordnet"); }
TEST(TestDPISOFAO, yeast) { filterTest("dpiso", "yeast"); }
TEST(TestDPISOFAO, youtube) { filterTest("dpiso", "youtube"); }

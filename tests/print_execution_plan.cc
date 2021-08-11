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
#include "gtest/gtest.h"

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/expand_edge_operator.h"
#include "ops/filters.h"
#include "ops/scans.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionConfig;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::Scan;
using circinus::NaivePlanner;
using circinus::NLFFilter;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::VertexID;

class PrintExecutionPlan : public testing::Test {
 protected:
  const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                              "patents", "wordnet", "yeast", "youtube"};

  std::vector<std::string> data_graph_paths_;
  std::vector<std::vector<std::string>> query_graph_paths_;

  void SetUp() override {
    data_graph_paths_.resize(datasets_.size());
    query_graph_paths_.reserve(datasets_.size());
    for (uint32_t dataset_i = 0; dataset_i < datasets_.size(); ++dataset_i) {
      auto& dataset = datasets_[dataset_i];
      data_graph_paths_[dataset_i] = dataset + "/data_graph/" + dataset + ".graph";
      query_graph_paths_.emplace_back();
      auto& paths = query_graph_paths_.back();
      uint32_t max_query_size = (dataset == "human" || dataset == "wordnet") ? 16 : 32;
      for (uint32_t i = 8; i <= max_query_size; i += 8) {
        paths.push_back(dataset + "/query_graph/query_dense_" + std::to_string(i) + "_1.graph");
        paths.push_back(dataset + "/query_graph/query_sparse_" + std::to_string(i) + "_1.graph");
      }
    }
  }

  std::vector<std::vector<VertexID>> getCandidateSets(const Graph& g, const QueryGraph& q) {
    std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
    ExecutionConfig config;
    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      config.setInputSize(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
      auto scan = circinus::Scan::newLDFScan(q.getVertexLabel(v), q.getVertexOutDegree(v), 0, config, 1);
      scan->addFilter(std::make_unique<NLFFilter>(&q, v));
      auto scan_ctx = scan->initScanContext(0);
      scan->scan(&g, &scan_ctx);
      candidates[v] = std::move(scan_ctx.candidates);
    }
    for (auto& candidate : candidates) {
      LOG(INFO) << "candidate set size " << candidate.size();
    }
    return candidates;
  }

  template <bool dynamic = false>
  void printExecutionPlan(uint32_t i) {
    Graph g(FLAGS_data_dir + "/" + data_graph_paths_[i]);  // load data graph
    for (uint32_t j = 0; j < query_graph_paths_[i].size(); ++j) {
      LOG(INFO) << "========================";
      LOG(INFO) << "graph " << data_graph_paths_[i] << " query " << query_graph_paths_[i][j];
      QueryGraph q(FLAGS_data_dir + "/" + query_graph_paths_[i][j]);  // load query graph
      auto candidates = getCandidateSets(g, q);                       // get candidates for each query vertex
      std::vector<double> candidate_cardinality;
      candidate_cardinality.reserve(candidates.size());
      for (auto& set : candidates) {
        candidate_cardinality.push_back(set.size());
      }
      NaivePlanner planner(&q, candidate_cardinality);
      ExecutionPlan* plan;
      if (dynamic) {
        plan = planner.generatePlanWithEagerDynamicCover();
      } else {
        plan = planner.generatePlan();
      }
      ASSERT_TRUE(plan != nullptr);
      plan->printPhysicalPlan();
    }
  }
};

TEST_F(PrintExecutionPlan, DBLP) { printExecutionPlan(0); }
TEST_F(PrintExecutionPlan, EU2005) { printExecutionPlan(1); }
TEST_F(PrintExecutionPlan, HPRD) { printExecutionPlan(2); }
TEST_F(PrintExecutionPlan, HUMAN) { printExecutionPlan(3); }
TEST_F(PrintExecutionPlan, PATENTS) { printExecutionPlan(4); }
TEST_F(PrintExecutionPlan, WORDNET) { printExecutionPlan(5); }
TEST_F(PrintExecutionPlan, YEAST) { printExecutionPlan(6); }
TEST_F(PrintExecutionPlan, YOUTUBE) { printExecutionPlan(7); }

TEST_F(PrintExecutionPlan, DynamicDBLP) { printExecutionPlan<true>(0); }
TEST_F(PrintExecutionPlan, DynamicEU2005) { printExecutionPlan<true>(1); }
TEST_F(PrintExecutionPlan, DynamicHPRD) { printExecutionPlan<true>(2); }
TEST_F(PrintExecutionPlan, DynamicHUMAN) { printExecutionPlan<true>(3); }
TEST_F(PrintExecutionPlan, DynamicPATENTS) { printExecutionPlan<true>(4); }
TEST_F(PrintExecutionPlan, DynamicWORDNET) { printExecutionPlan<true>(5); }
TEST_F(PrintExecutionPlan, DynamicYEAST) { printExecutionPlan<true>(6); }
TEST_F(PrintExecutionPlan, DynamicYOUTUBE) { printExecutionPlan<true>(7); }

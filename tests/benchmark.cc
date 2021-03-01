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

class Benchmark : public testing::Test {
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
    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      LDFScan scan(&q, v, &g);
      NLFFilter filter(&q, v);
      std::vector<VertexID> buffer;
      while (scan.Scan(&buffer, BATCH_SIZE) > 0) {
        filter.Filter(g, buffer, &candidates[v]);
        buffer.clear();
      }
    }
    for (auto& candidate : candidates) {
      LOG(INFO) << "candidate set size " << candidate.size();
    }
    return candidates;
  }

  void Execute(const Graph* g, const ExecutionPlan* plan) {
    LOG(INFO) << FLAGS_num_cores << " threads";
    ThreadPool threads(FLAGS_num_cores, plan);
    auto seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    if (plan->isInCover(plan->getRootQueryVertexID())) {
      for (size_t i = 0; i < seeds.size(); i += BATCH_SIZE) {
        size_t end = std::min(i + BATCH_SIZE, seeds.size());
        threads.addInitTask(0, std::vector<CompressedSubgraphs>(seeds.begin() + i, seeds.begin() + end), g);
      }
    } else {  // TODO(tatiana): this branch should not be reached
      threads.addInitTask(
          0, std::vector<CompressedSubgraphs>{CompressedSubgraphs(std::make_shared<std::vector<VertexID>>(seeds))}, g);
    }
    threads.start();
  }

  void run(uint32_t i) {
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
      NaivePlanner planner(&q, &candidate_cardinality);
      auto plan = planner.generatePlan();
      plan->setCandidateSets(candidates);  // swap
      plan->printPhysicalPlan();
      Execute(&g, plan);
    }
  }
};

TEST_F(Benchmark, DBLP) { run(0); }
TEST_F(Benchmark, EU2005) { run(1); }
TEST_F(Benchmark, HPRD) { run(2); }
TEST_F(Benchmark, HUMAN) { run(3); }
TEST_F(Benchmark, PATENTS) { run(4); }
TEST_F(Benchmark, WORDNET) { run(5); }
TEST_F(Benchmark, YEAST) { run(6); }
TEST_F(Benchmark, YOUTUBE) { run(7); }

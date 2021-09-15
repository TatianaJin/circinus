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

#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/filters.h"
#include "ops/scans.h"

using circinus::DUMMY_QUERY_VERTEX;
using circinus::ExecutionConfig;
using circinus::Graph;
using circinus::Scan;
using circinus::NLFFilter;
using circinus::QueryGraph;
using circinus::VertexID;
using circinus::WeightedBnB;

class TestFilterAndMWVC : public testing::Test {
 protected:
  const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                              "patents", "wordnet", "yeast", "youtube"};
  const std::string answer_path = "resources/test_filter_answer.txt";

  std::vector<std::string> data_graph_paths_;
  std::vector<std::vector<std::string>> query_graph_paths_;
  std::vector<std::vector<std::pair<std::vector<uint64_t>, std::vector<uint64_t>>>>
      candidate_size_fact_;  // { dataset : { query : { ldf, nlf } } }

  void SetUp() override {
    LOG(INFO) << "set up";
    std::ifstream input;
    input.open(answer_path);
    CHECK(input.is_open()) << "cannot open path " << answer_path << ". did you try cmake?";
    data_graph_paths_.resize(datasets_.size());
    query_graph_paths_.reserve(datasets_.size());
    std::string setting, filter;
    uint32_t num;
    candidate_size_fact_.resize(datasets_.size());
    for (uint32_t dataset_i = 0; dataset_i < datasets_.size(); ++dataset_i) {
      auto& dataset = datasets_[dataset_i];
      data_graph_paths_[dataset_i] = dataset + "/data_graph/" + dataset + ".graph";
      query_graph_paths_.emplace_back();
      auto& paths = query_graph_paths_.back();
      if (dataset == "human" || dataset == "wordnet") {
        paths.reserve(4);
        candidate_size_fact_[dataset_i].resize(4);
      } else {
        paths.reserve(8);
        candidate_size_fact_[dataset_i].resize(8);
      }
      uint32_t max_query_size = (dataset == "human" || dataset == "wordnet") ? 16 : 32;
      for (uint32_t i = 8; i <= max_query_size; i += 8) {
        paths.push_back(dataset + "/query_graph/query_dense_" + std::to_string(i) + "_1.graph");
        paths.push_back(dataset + "/query_graph/query_sparse_" + std::to_string(i) + "_1.graph");

        input >> setting >> filter;
        CHECK_EQ(setting, dataset + std::to_string(i) + "d1");
        CHECK_EQ(filter, "LDF");
        for (uint32_t c = 0; c < i; ++c) {
          input >> num;
          candidate_size_fact_[dataset_i][i / 4 - 2].first.push_back(num);
        }

        input >> setting >> filter;
        CHECK_EQ(setting, dataset + std::to_string(i) + "s1");
        CHECK_EQ(filter, "LDF");
        for (uint32_t c = 0; c < i; ++c) {
          input >> num;
          candidate_size_fact_[dataset_i][i / 4 - 1].first.push_back(num);
        }

        input >> setting >> filter;
        CHECK_EQ(setting, dataset + std::to_string(i) + "d1");
        CHECK_EQ(filter, "NLF");
        for (uint32_t c = 0; c < i; ++c) {
          input >> num;
          candidate_size_fact_[dataset_i][i / 4 - 2].second.push_back(num);
        }

        input >> setting >> filter;
        CHECK_EQ(setting, dataset + std::to_string(i) + "s1");
        CHECK_EQ(filter, "NLF");
        for (uint32_t c = 0; c < i; ++c) {
          input >> num;
          candidate_size_fact_[dataset_i][i / 4 - 1].second.push_back(num);
        }
      }
    }
    LOG(INFO) << "finished setup";
  }
};

TEST_F(TestFilterAndMWVC, FilterCorrectness) {
  for (uint32_t i = 0; i < datasets_.size(); ++i) {
    // load data graph
    Graph g(FLAGS_data_dir + "/" + data_graph_paths_[i]);
    LOG(INFO) << "load graph " << data_graph_paths_[i];
    for (uint32_t j = 0; j < query_graph_paths_[i].size(); ++j) {
      auto& q_path = query_graph_paths_[i][j];
      // load query graph
      QueryGraph q(FLAGS_data_dir + "/" + q_path);

      // get candidates for each query vertex
      auto start = std::chrono::high_resolution_clock::now();
      ExecutionConfig config;
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        config.setInputSize(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
        auto scan = Scan::newLDFScan(q.getVertexLabel(v), q.getVertexOutDegree(v), 0, config, 1);
        auto scan_ctx = scan->initScanContext(v, 0, {DUMMY_QUERY_VERTEX, 0});
        scan->scan(&g, &scan_ctx);
        size_t ldf_output_size = scan_ctx.candidates.size();
        std::vector<VertexID> candidates;
        NLFFilter filter(&q, v);
        filter.filter(g, scan_ctx.candidates, &candidates);
        EXPECT_EQ(candidate_size_fact_[i][j].first[v], ldf_output_size) << q_path;
        EXPECT_EQ(candidate_size_fact_[i][j].second[v], candidates.size()) << q_path;
      }
      auto end = std::chrono::high_resolution_clock::now();
      LOG(INFO) << q_path << ' ' << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
    }
  }
}

void PrintCover(const std::vector<int>& assignment) {
  std::stringstream ss;
  uint32_t size = 0;
  for (uint32_t i = 0; i < assignment.size(); ++i) {
    if (assignment[i] == 1) {
      ++size;
      ss << ' ' << i;
    }
  }
  LOG(INFO) << "vertex cover (size " << size << ')' << ss.str();
}

TEST_F(TestFilterAndMWVC, MWVC) {
  for (uint32_t i = 0; i < datasets_.size(); ++i) {
    // load data graph
    Graph g(FLAGS_data_dir + "/" + data_graph_paths_[i]);
    LOG(INFO) << "load graph " << data_graph_paths_[i];
    for (uint32_t j = 0; j < query_graph_paths_[i].size(); ++j) {
      auto& q_path = query_graph_paths_[i][j];
      // load query graph
      QueryGraph q(FLAGS_data_dir + "/" + q_path);
      // get candidates for each query vertex
      std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
      ExecutionConfig config;
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        config.setInputSize(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
        auto scan = Scan::newLDFScan(q.getVertexLabel(v), q.getVertexOutDegree(v), 0, config, 1);
        scan->addFilter(std::make_unique<NLFFilter>(&q, v));
        auto scan_ctx = scan->initScanContext(v, 0, {DUMMY_QUERY_VERTEX, 0});
        scan->scan(&g, &scan_ctx);
        candidates[v] = std::move(scan_ctx.candidates);
      }
      std::vector<double> vertex_weights;
      vertex_weights.reserve(candidates.size());
      for (auto& cs : candidates) {
        vertex_weights.push_back(cs.size());
      }
      WeightedBnB solver(&q, vertex_weights);
      solver.computeVertexCover();
      LOG(INFO) << solver.getBestCovers().size() << " best covers, weight" << solver.getBestObjective() << ' '
                << solver.getTimeToBest() << '/' << solver.getElapsedTime() << "ms";
      PrintCover(solver.getBestCovers().front());
    }
  }
}

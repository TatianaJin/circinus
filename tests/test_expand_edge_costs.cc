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

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/expand_edge_operator.h"
#include "ops/filters.h"
#include "ops/operators.h"
#include "ops/scans.h"
#include "utils/flags.h"
#include "utils/print_styles.h"

using circinus::CompressedSubgraphs;
using circinus::ExpandEdgeOperator;
using circinus::ExpandKeyToKeyVertexOperator;
using circinus::ExpandKeyToSetVertexOperator;
using circinus::ExpandSetToKeyVertexOperator;
using circinus::Graph;
using circinus::LDFScan;
using circinus::NLFFilter;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::VertexID;

#define BATCH_SIZE FLAGS_batch_size

DEFINE_string(data_dir, "/data/share/project/haxe/data/subgraph_matching_datasets", "The directory of datasets");

class TestExpandEdgeCosts : public testing::Test {
 protected:
  const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                              "patents", "wordnet", "yeast", "youtube"};

  std::vector<std::string> data_graph_paths_;
  std::vector<std::vector<std::string>> query_graph_paths_;

  void SetUp() override {
    LOG(INFO) << "BATCH_SIZE=" << BATCH_SIZE;
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

  static uint64_t getNumSubgraphs(const std::vector<CompressedSubgraphs>& outputs) {
    uint64_t n = 0;
    for (auto& output : outputs) {
      n += output.getNumSubgraphs();
    }
    return n;
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

  uint32_t getNumMatches(QueryVertexID parent, QueryVertexID target, const std::vector<int>& cover,
                         const std::vector<std::vector<VertexID>>& candidates,
                         const std::vector<CompressedSubgraphs>& seeds, const Graph& g) {
    std::unordered_map<QueryVertexID, uint32_t> indices;
    indices[parent] = 0;
    indices[target] = 1;
    auto op = ExpandEdgeOperator::newExpandEdgeOperator(parent, target, cover, indices);
    auto start = std::chrono::high_resolution_clock::now();
    op->setCandidateSets(&candidates[target]);
    op->input(seeds, &g);
    std::vector<CompressedSubgraphs> outputs;
    while (op->expand(&outputs, BATCH_SIZE) > 0) {
    }
    auto end = std::chrono::high_resolution_clock::now();
    delete op;
    auto ret = getNumSubgraphs(outputs);
    LOG(INFO) << ((cover[parent] == 1) ? "key" : "set") << "to" << ((cover[target] == 1) ? "key " : "set ") << parent
              << "->" << target << ' ' << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms " << ret << '/' << outputs.size();
    return ret;
  }

  uint32_t expandVertex(QueryVertexID parent, QueryVertexID target, const std::vector<int>& cover,
                        const std::vector<std::vector<VertexID>>& candidates,
                        const std::vector<CompressedSubgraphs>& seeds, const Graph& g) {
    std::unordered_map<QueryVertexID, uint32_t> indices;
    indices[parent] = 0;
    indices[target] = 1;

    // expand vertex operator
    circinus::TraverseOperator* op;
    if (cover[parent] == 1 && cover[target] == 1) {  // key to key
      op = new ExpandKeyToKeyVertexOperator(std::vector<QueryVertexID>{parent}, target, indices);
    } else if (cover[parent] == 1) {  // key to set
      op = new ExpandKeyToSetVertexOperator(std::vector<QueryVertexID>{parent}, target, indices);
    } else {  // set to key
      op = new ExpandSetToKeyVertexOperator(std::vector<QueryVertexID>{parent}, target, indices);
    }
    auto start = std::chrono::high_resolution_clock::now();
    op->setCandidateSets(&candidates[target]);
    op->input(seeds, &g);
    std::vector<CompressedSubgraphs> outputs;
    while (op->expand(&outputs, BATCH_SIZE) > 0) {
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ret = getNumSubgraphs(outputs);
    outputs.clear();
    auto expand_vertex_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    delete op;

    // expand vertex operator
    auto op_expand_edge = ExpandEdgeOperator::newExpandEdgeOperator(parent, target, cover, indices);
    start = std::chrono::high_resolution_clock::now();
    op_expand_edge->setCandidateSets(&candidates[target]);
    op_expand_edge->input(seeds, &g);
    while (op_expand_edge->expand(&outputs, BATCH_SIZE) > 0) {
    }
    end = std::chrono::high_resolution_clock::now();
    EXPECT_EQ(ret, getNumSubgraphs(outputs));
    auto expand_edge_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    delete op_expand_edge;

    if (expand_vertex_time < expand_edge_time) {
      LOG(INFO) << ((cover[parent] == 1) ? "key" : "set") << "to" << ((cover[target] == 1) ? "key " : "set ") << parent
                << "->" << target << " vertex/edge " << GREEN(expand_vertex_time << '/' << expand_edge_time) << "ms "
                << ret << '/' << outputs.size();
    } else if (expand_vertex_time > expand_edge_time) {
      LOG(INFO) << ((cover[parent] == 1) ? "key" : "set") << "to" << ((cover[target] == 1) ? "key " : "set ") << parent
                << "->" << target << " vertex/edge " << RED(expand_vertex_time << '/' << expand_edge_time) << "ms "
                << ret << '/' << outputs.size();
    } else {
      LOG(INFO) << ((cover[parent] == 1) ? "key" : "set") << "to" << ((cover[target] == 1) ? "key " : "set ") << parent
                << "->" << target << " vertex/edge " << expand_vertex_time << '/' << expand_edge_time << "ms " << ret
                << '/' << outputs.size();
    }
    return ret;
  }

  template <bool compare_expand_vertex = false>
  void expandEdges(uint32_t i) {
    Graph g(FLAGS_data_dir + "/" + data_graph_paths_[i]);  // load data graph
    for (uint32_t j = 0; j < query_graph_paths_[i].size(); ++j) {
      LOG(INFO) << "========================";
      LOG(INFO) << "graph " << data_graph_paths_[i] << " query " << query_graph_paths_[i][j];
      QueryGraph q(FLAGS_data_dir + "/" + query_graph_paths_[i][j]);  // load query graph
      auto candidates = getCandidateSets(g, q);                       // get candidates for each query vertex
      // expand the first edge from each v in q in different ways
      std::vector<int> cover(q.getNumVertices(), 0);
      std::vector<CompressedSubgraphs> seed_keys, seed_sets;
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        // generate seed graphs
        seed_keys.clear();
        seed_keys.reserve(candidates[v].size());
        for (auto c : candidates[v]) {
          seed_keys.emplace_back(c);
        }
        seed_sets.clear();
        seed_sets.emplace_back(std::make_shared<std::vector<VertexID>>(candidates[v]));
        // expand from v to each neighbor of v
        auto neighbors = q.getOutNeighbors(v);
        for (uint32_t nb = 0; nb < neighbors.second; ++nb) {
          auto target = neighbors.first[nb];
          uint32_t n_output;
          { /* key to set */
            cover[v] = 1;
            cover[target] = 0;
            if (compare_expand_vertex) {
              n_output = expandVertex(v, target, cover, candidates, seed_keys, g);
            } else {
              n_output = getNumMatches(v, target, cover, candidates, seed_keys, g);
            }
          }
          { /* key to key */
            cover[v] = 1;
            cover[target] = 1;
            if (compare_expand_vertex) {
              EXPECT_EQ(n_output, expandVertex(v, target, cover, candidates, seed_keys, g));
            } else {
              EXPECT_EQ(n_output, getNumMatches(v, target, cover, candidates, seed_keys, g));
            }
          }
          { /* set to key */
            cover[v] = 0;
            cover[target] = 1;
            if (compare_expand_vertex) {
              ASSERT_EQ(n_output, expandVertex(v, target, cover, candidates, seed_sets, g));
            } else {
              ASSERT_EQ(n_output, getNumMatches(v, target, cover, candidates, seed_sets, g));
            }
          }
          LOG(INFO) << "-------";
        }
      }
    }
  }
};

TEST_F(TestExpandEdgeCosts, DBLP) { expandEdges(0); }
TEST_F(TestExpandEdgeCosts, EU2005) { expandEdges(1); }
TEST_F(TestExpandEdgeCosts, HPRD) { expandEdges(2); }
TEST_F(TestExpandEdgeCosts, HUMAN) { expandEdges(3); }
TEST_F(TestExpandEdgeCosts, PATENTS) { expandEdges(4); }
TEST_F(TestExpandEdgeCosts, WORDNET) { expandEdges(5); }
TEST_F(TestExpandEdgeCosts, YEAST) { expandEdges(6); }
TEST_F(TestExpandEdgeCosts, YOUTUBE) { expandEdges(7); }

TEST_F(TestExpandEdgeCosts, CompareDBLP) { expandEdges<true>(0); }
TEST_F(TestExpandEdgeCosts, CompareEU2005) { expandEdges<true>(1); }
TEST_F(TestExpandEdgeCosts, CompareHPRD) { expandEdges<true>(2); }
TEST_F(TestExpandEdgeCosts, CompareHUMAN) { expandEdges<true>(3); }
TEST_F(TestExpandEdgeCosts, ComparePATENTS) { expandEdges<true>(4); }
TEST_F(TestExpandEdgeCosts, CompareWORDNET) { expandEdges<true>(5); }
TEST_F(TestExpandEdgeCosts, CompareYEAST) { expandEdges<true>(6); }
TEST_F(TestExpandEdgeCosts, CompareYOUTUBE) { expandEdges<true>(7); }

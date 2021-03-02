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
#include <numeric>
#include <string>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/graph.h"
#include "ops/expand_edge_operator.h"
#include "ops/expand_into_operator.h"
#include "ops/expand_key_to_key_vertex_operator.h"
#include "ops/expand_key_to_set_vertex_operator.h"
#include "ops/expand_set_to_key_vertex_operator.h"
#include "ops/filters.h"
#include "ops/operator.h"
#include "ops/scans.h"
#include "plan/execution_plan.h"

using circinus::CompressedSubgraphs;
using circinus::ExpandKeyToSetVertexOperator;
using circinus::ExpandSetToKeyVertexOperator;
using circinus::ExpandKeyToKeyVertexOperator;
using circinus::ExpandEdgeOperator;
using circinus::ExpandIntoOperator;
using circinus::Graph;
using circinus::LDFScan;
using circinus::NLFFilter;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::VertexID;
using circinus::WeightedBnB;
using circinus::Operator;

#define BATCH_SIZE 32

DEFINE_string(data_dir, "/data/share/project/haxe/data/subgraph_matching_datasets", "The directory of datasets");

class TestSubgraphMatching : public testing::Test {
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

  uint64_t getNumSubgraphs(const std::vector<CompressedSubgraphs>& outputs) {
    uint64_t n = 0;
    for (auto& output : outputs) {
      n += output.getNumSubgraphs();
    }
    return n;
  }

  void getMatchingOrder(std::vector<QueryVertexID>& matching_order, std::vector<bool>& visited, QueryGraph* q,
                        QueryVertexID uid) {
    if (visited[uid]) {
      return;
    }
    visited[uid] = true;
    matching_order.push_back(uid);
    const auto& neighbors = q->getOutNeighbors(uid);
    for (uint32_t i = 0; i < neighbors.second; ++i) {
      auto vid = neighbors.first[i];
      getMatchingOrder(matching_order, visited, q, vid);
    }
  }

  void outputLog(const std::string name, std::vector<QueryVertexID> parents, QueryVertexID target, uint32_t time_usage,
                 uint32_t output_num) {
    std::string parents_string;
    for (auto uid : parents) {
      parents_string += std::to_string(uid) + " ";
    }
    LOG(INFO) << name << " " << parents_string << "->" << target << ' ' << time_usage << "ms " << output_num;
  }

  void subgraphMatching(uint32_t graph_id) {
    // load data graph
    Graph g(FLAGS_data_dir + "/" + data_graph_paths_[graph_id]);
    LOG(INFO) << "load graph " << data_graph_paths_[graph_id];
    for (uint32_t query_id = 0; query_id < query_graph_paths_[graph_id].size(); ++query_id) {
      auto& q_path = query_graph_paths_[graph_id][query_id];
      // load query graph
      QueryGraph q(FLAGS_data_dir + "/" + q_path);
      LOG(INFO) << "query graph " << q_path;
      // get candidates for each query vertex
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

      std::vector<double> vertex_weights;
      vertex_weights.reserve(candidates.size());
      for (auto& cs : candidates) {
        vertex_weights.push_back(cs.size());
      }
      WeightedBnB solver(&q, vertex_weights);
      solver.computeVertexCover();
      LOG(INFO) << solver.getBestCovers().size() << " best covers, weight" << solver.getBestObjective() << ' '
                << solver.getTimeToBest() << '/' << solver.getElapsedTime() << "s";
      auto cover = solver.getBestCovers().front();
      for (auto& set : candidates) {
        LOG(INFO) << "candidate set size " << set.size();
      }

      std::unordered_map<QueryVertexID, uint32_t> indices;

      std::vector<QueryVertexID> matching_order;
      std::vector<bool> visited(q.getNumVertices(), false);
      getMatchingOrder(matching_order, visited, &q, 0);

      for (uint32_t i = 0; i < matching_order.size(); i++) {
        LOG(INFO) << "qid " << matching_order[i] << " cover " << cover[matching_order[i]];
        cover[matching_order[i]] = cover[matching_order[i]] > 0;
      }

      indices[matching_order[0]] = 0;
      uint32_t n_keys = cover[matching_order[0]];
      uint32_t n_sets = 1 - n_keys, batch_size = 10000;
      std::unordered_set<QueryVertexID> existing_vertices;
      existing_vertices.insert(matching_order[0]);
      std::vector<std::vector<CompressedSubgraphs>> outputs;
      outputs.resize(2);
      if (cover[matching_order[0]]) {
        for (auto c : candidates[matching_order[0]]) {
          outputs[0].emplace_back(c);
        }
      } else {
        outputs[0].emplace_back(std::make_shared<std::vector<VertexID>>(std::move(candidates[matching_order[0]])));
      }

      std::vector<QueryVertexID> key_parents;
      std::vector<QueryVertexID> set_parents;
      for (uint32_t i = 0, j = 1; j < matching_order.size(); ++i, ++j) {
        uint32_t target = matching_order[j];
        outputs[!(i & 1)].clear();

        // record query vertex index
        if (cover[target] == 1) {
          indices[target] = n_keys;
          ++n_keys;
        } else {
          indices[target] = n_sets;
          ++n_sets;
        }

        auto neighbors = q.getOutNeighbors(target);
        for (uint32_t j = 0; j < neighbors.second; ++j) {
          if (existing_vertices.count(neighbors.first[j]) != 0) {
            if (cover[neighbors.first[j]]) {
              key_parents.push_back(neighbors.first[j]);
            } else {
              set_parents.push_back(neighbors.first[j]);
            }
          }
        }
        if (key_parents.size() + set_parents.size() == 1) {  // only one parent, ExpandEdge
          auto front = key_parents.size() == 1 ? key_parents.front() : set_parents.front();
          auto start = std::chrono::high_resolution_clock::now();
          auto op = ExpandEdgeOperator::newExpandEdgeOperator(front, target, cover, indices);
          op->setCandidateSets(&candidates[target]);
          op->input(outputs[i & 1], &g);
          while (op->expand(&outputs[!(i & 1)], batch_size) > 0) {
          }
          auto end = std::chrono::high_resolution_clock::now();
          outputLog("expand edge", key_parents.size() == 1 ? key_parents : set_parents, target,
                    std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
                    outputs[!(i & 1)].size());
          delete op;
        } else {  // more than one parents, ExpandVertex (use set intersection)
          if (cover[target]) {
            if (key_parents.size() != 0 && set_parents.size() != 0) {
              // key to key
              LOG(INFO) << "key and set to key";
              auto start = std::chrono::high_resolution_clock::now();
              auto op0 = new ExpandKeyToKeyVertexOperator(key_parents, target, indices);
              op0->setCandidateSets(&candidates[target]);
              op0->input(outputs[i & 1], &g);
              while (op0->expand(&outputs[!(i & 1)], batch_size) > 0) {
              }
              auto end = std::chrono::high_resolution_clock::now();
              outputLog("key to key", key_parents, target,
                        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
                        outputs[!(i & 1)].size());
              i++;
              outputs[!(i & 1)].clear();
              delete op0;
              start = std::chrono::high_resolution_clock::now();
              auto op1 = new ExpandIntoOperator(set_parents, target, indices);
              op1->setCandidateSets(&candidates[target]);
              op1->input(outputs[i & 1], &g);
              while (op1->expand(&outputs[!(i & 1)], batch_size) > 0) {
              }
              end = std::chrono::high_resolution_clock::now();
              outputLog("expand into", set_parents, target,
                        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
                        outputs[!(i & 1)].size());
              delete op1;
            } else if (key_parents.size() != 0) {
              auto start = std::chrono::high_resolution_clock::now();
              auto op = new ExpandKeyToKeyVertexOperator(key_parents, target, indices);
              op->setCandidateSets(&candidates[target]);
              op->input(outputs[i & 1], &g);
              while (op->expand(&outputs[!(i & 1)], batch_size) > 0) {
              }
              auto end = std::chrono::high_resolution_clock::now();
              outputLog("key to key", key_parents, target,
                        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
                        outputs[!(i & 1)].size());
              delete op;
            } else {
              LOG(INFO) << "set to key";
              auto start = std::chrono::high_resolution_clock::now();
              auto op = new ExpandSetToKeyVertexOperator(set_parents, target, indices);
              for (auto uid : set_parents) {
                LOG(INFO) << uid;
              }
              op->setCandidateSets(&candidates[target]);
              op->input(outputs[i & 1], &g);
              LOG(INFO) << target;
              while (op->expand(&outputs[!(i & 1)], batch_size) > 0) {
              }
              auto end = std::chrono::high_resolution_clock::now();
              outputLog("set to key", set_parents, target,
                        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
                        outputs[!(i & 1)].size());
              delete op;
            }
          } else {
            LOG(INFO) << "key to set";
            auto start = std::chrono::high_resolution_clock::now();
            auto op = new ExpandKeyToSetVertexOperator(key_parents, target, indices);
            op->setCandidateSets(&candidates[target]);
            op->input(outputs[i & 1], &g);
            while (op->expand(&outputs[!(i & 1)], batch_size) > 0) {
            }
            auto end = std::chrono::high_resolution_clock::now();
            outputLog("key to set", key_parents, target,
                      std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(),
                      outputs[!(i & 1)].size());
            delete op;
          }
        }
        LOG(INFO) << " matching number " << getNumSubgraphs(outputs[!(i & 1)]);
        key_parents.clear();
        set_parents.clear();
        existing_vertices.insert(target);
      }
    }
  }
};

TEST_F(TestSubgraphMatching, DBLP) { subgraphMatching(0); }
TEST_F(TestSubgraphMatching, EU2005) { subgraphMatching(1); }
TEST_F(TestSubgraphMatching, HPRD) { subgraphMatching(2); }
TEST_F(TestSubgraphMatching, HUMAN) { subgraphMatching(3); }
TEST_F(TestSubgraphMatching, PATENTS) { subgraphMatching(4); }
TEST_F(TestSubgraphMatching, WORDNET) { subgraphMatching(5); }
TEST_F(TestSubgraphMatching, YEAST) { subgraphMatching(6); }
TEST_F(TestSubgraphMatching, YOUTUBE) { subgraphMatching(7); }

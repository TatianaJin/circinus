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

#include <fstream>
#include <string>

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "graph/query_graph.h"
#include "ops/filters.h"

using circinus::VertexID;
using circinus::Graph;

class TestGraph : public testing::Test {};

TEST_F(TestGraph, Load) {
  Graph g("resources/human.graph");
  EXPECT_EQ(g.getNumVertices(), 4674u);
  EXPECT_EQ(g.getNumEdges(), 86282u);
  EXPECT_EQ(g.getGraphMaxDegree(), 771u);
  auto labels = g.getLabels();
  EXPECT_EQ(labels.size(), 44u);
  uint64_t max_label_frequency = 0;
  for (auto l : labels) {
    max_label_frequency = std::max(max_label_frequency, g.getVertexCardinalityByLabel(l));
  }
  EXPECT_EQ(max_label_frequency, 683u);
}

TEST_F(TestGraph, GetVerticesByLabel) {
  Graph g("resources/human.graph");
  auto labels = g.getLabels();
  for (auto l : labels) {
    EXPECT_EQ(g.getVerticesByLabel(l)->size(), g.getVertexCardinalityByLabel(l));
  }
}

TEST_F(TestGraph, BinarySerDe) {
  std::string graph_path = "resources/human.graph";

  Graph g(graph_path);
  g.saveAsBinary(graph_path + ".bin");
  Graph b;
  std::ifstream input(graph_path + ".bin", std::ios::binary);
  b.loadCompressed(input);

  ASSERT_EQ(g.getNumVertices(), b.getNumVertices());
  ASSERT_EQ(g.getNumEdges(), b.getNumEdges());
  ASSERT_EQ(g.getGraphMaxDegree(), b.getGraphMaxDegree());

  auto g_labels = g.getLabels();
  auto b_labels = b.getLabels();
  ASSERT_EQ(g_labels.size(), b_labels.size());
  std::sort(g_labels.begin(), g_labels.end());
  std::sort(b_labels.begin(), b_labels.end());
  for (uint32_t i = 0; i < b_labels.size(); ++i) {
    ASSERT_EQ(g_labels[i], b_labels[i]);
  }

  for (VertexID i = 0; i < g.getNumVertices(); ++i) {
    ASSERT_EQ(g.getVertexOutDegree(i), b.getVertexOutDegree(i));
    ASSERT_EQ(g.getVertexLabel(i), b.getVertexLabel(i));
    auto g_nbs = g.getOutNeighbors(i);
    auto b_nbs = b.getOutNeighbors(i);
    EXPECT_EQ(g_nbs.second, b_nbs.second);
    for (uint32_t j = 0; j < g_nbs.second; ++j) {
      ASSERT_EQ(g_nbs.first[j], b_nbs.first[j]);
    }
  }

  for (auto label : g.getLabels()) {
    ASSERT_EQ(g.getVertexCardinalityByLabel(label), b.getVertexCardinalityByLabel(label));
    auto* gv = g.getVerticesByLabel(label);
    auto* bv = g.getVerticesByLabel(label);
    ASSERT_EQ(gv->size(), bv->size());
    ASSERT_EQ(g.getNumVerticesByLabel(label), b.getNumVerticesByLabel(label));
    for (uint32_t i = 0; i < gv->size(); ++i) {
      ASSERT_EQ((*gv)[i], (*bv)[i]) << "label " << label << " idx " << i;
    }
  }
}

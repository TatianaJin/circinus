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
#include "graph/bipartite_graph.h"
#include "graph/query_graph.h"
#include "gtest/gtest.h"
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

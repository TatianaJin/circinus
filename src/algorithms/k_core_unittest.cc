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

#include "algorithms/k_core.h"

#include "gtest/gtest.h"

#include "graph/query_graph.h"

using circinus::QueryGraph;
using circinus::TwoCoreSolver;

class TestTwoCoreSolver : public testing::Test {};

TEST_F(TestTwoCoreSolver, Correctness) {
  QueryGraph graph("resources/test_vertex_cover1.graph");
  TwoCoreSolver solver(&graph);
  solver.get2CoreVertices();
  // TODO(tatiana): check answer
}

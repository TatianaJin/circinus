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
}

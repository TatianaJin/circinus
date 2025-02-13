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

#include <string>
#include <vector>

#include "glog/logging.h"
#include "gtest/gtest.h"

#include "algorithms/minimum_weight_vertex_cover.h"
#include "graph/query_graph.h"

using circinus::WeightedBnB;
using circinus::QueryGraph;
using circinus::QueryVertexID;

class TestMinimumWeightVertexCoverAlgorithm : public testing::Test {
 protected:
  std::vector<std::string> test_graph_path_ = {
      "resources/test_vertex_cover1.graph", "resources/test_vertex_cover2.graph", "resources/test_vertex_cover3.graph",
  };
  std::vector<std::pair<uint32_t, uint32_t>> test_graph_answers_ = {
      {17, 144}, {27, 819}, {4, 5}};  // {best_cover_size, num_best_covers}
};

TEST_F(TestMinimumWeightVertexCoverAlgorithm, InitAndDelete) {
  QueryGraph g;
  std::vector<double> vertex_weights(g.getNumVertices());
  auto solver = new WeightedBnB(&g, vertex_weights);
  EXPECT_TRUE(solver != nullptr);
  delete solver;
}

TEST_F(TestMinimumWeightVertexCoverAlgorithm, Correctness) {
  for (uint32_t i = 0; i < test_graph_path_.size(); ++i) {
    QueryGraph g(test_graph_path_[i]);
    std::vector<double> vertex_weights(g.getNumVertices(), 1);
    WeightedBnB solver(&g, vertex_weights, 120);  // the test case should finish within 2 minutes
    solver.computeVertexCover();
    LOG(INFO) << "Elapsed time " << solver.getElapsedTime() << " first result found at " << solver.getTimeToBest();
    LOG(INFO) << "minimum vertex cover weight " << solver.getBestObjective() << " number "
              << solver.getBestCovers().size();
    ASSERT_EQ(test_graph_answers_[i].first, solver.getBestObjective());
    ASSERT_EQ(test_graph_answers_[i].second, solver.getBestCovers().size());
    for (auto& assignment : solver.getBestCovers()) {
      ASSERT_EQ(assignment.size(), g.getNumVertices());

      // check the number of vertices assigned to be in the cover
      int vertex_cover_size = 0;
      std::stringstream ss;
      int idx = 0;
      for (auto flag : assignment) {
        vertex_cover_size += (flag == 1);
        if (flag == 1) ss << " " << idx;
        ++idx;
      }
      // DLOG(INFO) << "vertex cover" << ss.str();
      EXPECT_EQ(vertex_cover_size, solver.getBestObjective());

      // check if it is a vertex cover
      for (QueryVertexID i = 0; i < g.getNumVertices(); ++i) {
        if (assignment[i]) {  // vertex i is in the cover, all its out edges are covered
          continue;
        }
        auto neighbors = g.getOutNeighbors(i);
        // the destination vertex should be in the cover to cover j-th out-edge of i
        for (uint32_t j = 0; j < neighbors.second; ++j) {
          EXPECT_EQ(1, assignment[neighbors.first[j]]);
        }
      }
    }
  }
}

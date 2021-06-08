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

#include <algorithm>
#include <fstream>
#include <string>

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "graph/graph.h"
#include "graph/partitioned_graph.h"
#include "ops/filters.h"

using circinus::Graph;
using circinus::ReorderedPartitionedGraph;
using circinus::VertexID;

class TestReorderedPartitionedGraph : public testing::Test {
 protected:
  uint32_t n_partitions_ = 10;

  template <class GraphA, class GraphB>
  void checkLabel(const GraphA& g, const GraphB& b) {
    auto g_labels = g.getLabels();
    auto b_labels = b.getLabels();
    ASSERT_EQ(g_labels.size(), b_labels.size());
    std::sort(g_labels.begin(), g_labels.end());
    std::sort(b_labels.begin(), b_labels.end());
    for (uint32_t i = 0; i < b_labels.size(); ++i) {
      ASSERT_EQ(g_labels[i], b_labels[i]);
    }
  }

  void checkGraphConsistency(const ReorderedPartitionedGraph& g, const Graph& b) {
    ASSERT_EQ(g.getNumVertices(), b.getNumVertices());
    ASSERT_EQ(g.getNumEdges(), b.getNumEdges());
    ASSERT_EQ(g.getGraphMaxDegree(), b.getGraphMaxDegree());

    checkLabel(g, b);

    for (VertexID i = 0; i < g.getNumVertices(); ++i) {
      ASSERT_EQ(g.getVertexOutDegree(i), b.getVertexOutDegree(g.getOriginalVertexId(i)));
      ASSERT_EQ(g.getVertexLabel(i), b.getVertexLabel(g.getOriginalVertexId(i)));
      auto g_nbs = g.getOutNeighbors(i);
      auto b_nbs = b.getOutNeighbors(g.getOriginalVertexId(i));
      EXPECT_EQ(g_nbs.second, b_nbs.second);
      std::vector<VertexID> g_nbrs_sorted_original_id(g_nbs.second);
      for (uint32_t j = 0; j < g_nbs.second; ++j) {
        g_nbrs_sorted_original_id[j] = g.getOriginalVertexId(g_nbs.first[j]);
      }
      std::sort(g_nbrs_sorted_original_id.begin(), g_nbrs_sorted_original_id.end());
      for (uint32_t j = 0; j < g_nbs.second; ++j) {
        ASSERT_EQ(g_nbrs_sorted_original_id[j], b_nbs.first[j]);
      }
    }

    std::vector<bool> label_mask(g.getNumVertices(), false);
    for (auto label : g.getLabels()) {
      ASSERT_EQ(g.getVertexCardinalityByLabel(label), b.getVertexCardinalityByLabel(label));
      VertexID label_frequency = 0;
      for (uint32_t i = 0; i < n_partitions_; ++i) {
        auto gv = g.getVertexRangeByLabel(label, i);
        label_frequency += gv.second - gv.first;
        for (auto j = gv.first; j < gv.second; ++j) {
          label_mask[g.getOriginalVertexId(j)] = true;
        }
      }
      auto bv = b.getVerticesByLabel(label);
      ASSERT_EQ(label_frequency, bv->size());
      ASSERT_EQ(g.getVertexCardinalityByLabel(label), b.getVertexCardinalityByLabel(label));
      for (uint32_t i = 0; i < label_frequency; ++i) {
        ASSERT_TRUE(label_mask[(*bv)[i]]) << "label " << label << " idx " << i;
      }
    }
  }
};

TEST_F(TestReorderedPartitionedGraph, ConsistencyWithGraph) {
  std::string graph_path = "resources/human.graph";
  ReorderedPartitionedGraph g(graph_path, n_partitions_);
  Graph b(graph_path);
  for (uint32_t partition = 0; partition < n_partitions_; ++partition) {
    auto range = g.getPartitionRange(partition);
    for (auto v = range.first; v < range.second; ++v) {
      auto glabel = g.getVertexLabelInPartition(v, partition);
      auto blabel = b.getVertexLabel(g.getOriginalVertexId(v));
      ASSERT_EQ(glabel, blabel) << "g.v=" << v << " original=" << g.getOriginalVertexId(v);
    }
  }
}

TEST_F(TestReorderedPartitionedGraph, SerdeConsistencyWithGraph) {
  std::string graph_path = "resources/human.graph";
  ReorderedPartitionedGraph g(graph_path, n_partitions_);
  g.dumpToFile(graph_path + ".check");
  Graph b(graph_path + ".check");
  checkGraphConsistency(g, b);
}

TEST_F(TestReorderedPartitionedGraph, ConversionFromGraph) {
  std::string graph_path = "resources/human.graph";
  Graph b(graph_path);
  ReorderedPartitionedGraph g(b, n_partitions_);
  checkGraphConsistency(g, b);
}

TEST_F(TestReorderedPartitionedGraph, BinarySerDe) {
  std::string graph_path = "resources/human.graph";

  ReorderedPartitionedGraph g(graph_path);
  g.saveAsBinary(graph_path + ".bin");
  ReorderedPartitionedGraph b;
  std::ifstream input(graph_path + ".bin", std::ios::binary);
  b.loadUndirectedGraphFromBinary(input);

  ASSERT_EQ(g.getNumVertices(), b.getNumVertices());
  ASSERT_EQ(g.getNumEdges(), b.getNumEdges());
  ASSERT_EQ(g.getGraphMaxDegree(), b.getGraphMaxDegree());

  checkLabel(g, b);

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
    for (uint32_t i = 0; i < n_partitions_; ++i) {
      auto gv = g.getVertexRangeByLabel(label, i);
      auto bv = b.getVertexRangeByLabel(label, i);
      ASSERT_EQ(gv.second - gv.first, bv.second - bv.first);
      ASSERT_EQ(g.getVertexCardinalityByLabel(label), b.getVertexCardinalityByLabel(label));
      auto frequency = gv.second - gv.first;
      for (uint32_t i = 0; i < frequency; ++i) {
        ASSERT_EQ(gv.first + i, bv.first + i) << "label " << label << " idx " << i;
      }
    }
  }
}

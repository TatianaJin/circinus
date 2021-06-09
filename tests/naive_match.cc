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

#include "gflags/gflags.h"
#include "glog/logging.h"
#ifdef WITH_GPERF
#include "gperftools/profiler.h"
#endif
#include "gtest/gtest.h"

#include "exec/thread_pool.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters.h"
#include "ops/operators.h"
#include "ops/order.h"
#include "ops/scans.h"
#include "ops/types.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/hashmap.h"
#include "utils/profiler.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionConfig;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::GraphType;
using circinus::NaivePlanner;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::CoverNode;
using circinus::QueryType;
using circinus::TraverseOperator;
using circinus::VertexID;


DEFINE_string(naive_datagraph, "/data/share/users/qlma/query-graph-output/query_dense_4_1.graph",
              "data graph file path");
DEFINE_string(naive_querygraph, "/data/share/users/qlma/query-graph-output/query_dense_4_1.graph",
              "query graph file path");


  void getDFSOrder(std::vector<QueryVertexID>& use_order, QueryGraph q, QueryVertexID i, std::vector<bool>& visited) {
    use_order.push_back(i);
    visited[i] = 1;
    auto[vec, len] = q.getOutNeighbors(i);
    for (uint32_t j = 0; j < len; ++j) {
      QueryVertexID nextv = vec[j];
      if (!visited[nextv]) getDFSOrder(use_order, q, nextv, visited);
    }
  }
  std::vector<std::vector<VertexID>> naiveGetCandidateSets(const Graph& g, const QueryGraph& q) {
    std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
    ExecutionConfig config;
    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      config.setInputSize(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
      if (config.getInputSize() <= 0) break;
      auto scan = circinus::Scan::newLDFScan(q.getVertexLabel(v), q.getVertexOutDegree(v), 0, config, 1);
      auto scan_ctx = scan->initScanContext(0);
      scan->scan(&g, &scan_ctx);
      candidates[v] = std::move(scan_ctx.candidates);
      // LOG(INFO) << "query vertex " << v << ' ' << scan->toString();
    }
    return candidates;
  }
  void naiverun(const std::string& datagraph, const std::string& naive_querygraph) {
    Graph g(datagraph);
    QueryGraph q(naive_querygraph);
    if(g.getNumEdges()!=q.getNumEdges())return;
    if(g.getNumVertices()!=q.getNumVertices())return;
    std::vector<QueryVertexID> use_order;
    std::vector<bool> visited(q.getNumVertices());
    getDFSOrder(use_order, q, 0, visited);
    auto candidates = naiveGetCandidateSets(g, q);
    std::vector<double> candidate_cardinality;
    candidate_cardinality.reserve(candidates.size());
    for (auto& set : candidates) {
      double size = set.size();
      candidate_cardinality.push_back(size);
    }
    NaivePlanner planner(&q, &candidate_cardinality, GraphType::Normal);
    ExecutionPlan* plan = planner.generatePlanWithoutCompression(use_order);
    plan->setCandidateSets(candidates);
    plan->printPhysicalPlan();
    plan->getOutputs().init(1).limit(1);
		auto seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    if (plan->isInCover(plan->getRootQueryVertexID()))
    {
        plan->getOperators().execute(&g, std::vector<CompressedSubgraphs>(seeds.begin(), seeds.end()), 0);
    }
    else
    {
        std::vector<CompressedSubgraphs> input;
      input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(seeds)));
      plan->getOperators().execute(&g, input, 0);
    }
    auto n_matches = plan->getOutputs().getCount();
    if (n_matches) std::cout << "MATCH!!!\n";
  }

int main(int argc, char** argv) {

    #ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_log_dir == "") {
    google::LogToStderr();
  }

naiverun(FLAGS_naive_datagraph, FLAGS_naive_querygraph);
}

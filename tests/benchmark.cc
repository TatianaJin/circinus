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
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

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
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/profiler.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::LDFScan;
using circinus::NaivePlanner;
using circinus::NLFFilter;
using circinus::CFLFilter;
using circinus::CFLOrder;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::Task;
using circinus::ThreadPool;
using circinus::VertexID;
using circinus::Profiler;

#define BATCH_SIZE FLAGS_batch_size
#define toSeconds(start, end) \
  (((double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1e9)

DEFINE_string(data_dir, "/data/share/project/haxe/data/subgraph_matching_datasets", "The directory of datasets");
DEFINE_string(output_file, "", "The output file path");
DEFINE_string(dataset, "dblp", "The dataset to use");
DEFINE_string(query_mode, "dense", "Dense or sparse query");
DEFINE_uint64(query_size, 8, "The query size.");
DEFINE_uint64(match_limit, 1e5, "The limit of matches to find");
DEFINE_uint64(query_index, 1, "The index of query in the same category");
DEFINE_string(match_order, "", "Matching order");
DEFINE_string(filter, "nlf", "filter");
DEFINE_string(profile_file, "", "profile file");
DEFINE_string(vertex_cover, "static", "Vertex cover strategy: static, eager");

enum VertexCoverStrategy : uint32_t { Static = 0, Eager };

class Benchmark {
 protected:
  const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                              "patents", "wordnet", "yeast", "youtube"};

 public:
  template <VertexCoverStrategy vcs>
  void run(const std::string& dataset, uint32_t query_size, const std::string& query_mode, uint32_t index,
           std::ostream* out) {
    auto graph_path = dataset + "/data_graph/" + dataset + ".graph";
    auto query_path = dataset + "/query_graph/query_" + query_mode + "_" + std::to_string(query_size) + "_" +
                      std::to_string(index) + ".graph";

    (*out) << dataset << ',' << query_size << ',' << query_mode << ',' << index << ',';
    run<vcs>(graph_path, query_path, out);
  }

 protected:
  std::vector<std::vector<VertexID>> getCandidateSets(const Graph& g, const QueryGraph& q) {
    std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
    std::vector<uint32_t> candidate_size(q.getNumVertices());
    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      candidates[v].reserve(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
      LDFScan scan(&q, v, &g);
      NLFFilter filter(&q, v);
      std::vector<VertexID> buffer;
      buffer.reserve(BATCH_SIZE);
      while (scan.Scan(&buffer, BATCH_SIZE) > 0) {
        filter.Filter(g, buffer, &candidates[v]);
        buffer.clear();
      }
      candidate_size[v] = candidates[v].size();
    }
    if (FLAGS_filter == "cfl") {
      CFLOrder cfl_order;
      QueryVertexID start_vertex = cfl_order.getStartVertex(&g, &q, candidate_size);
      LOG(INFO) << "cfl order get start vertex " << start_vertex;
      CFLFilter cfl_filter(&q, &g, start_vertex);
      cfl_filter.Filter(candidates);
    }
    for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
      LOG(INFO) << "vertex " << v << " " << candidate_size[v] << "/" << candidates[v].size();
    }
    return candidates;
  }

  uint64_t getNumIsomorphicSubgraphs(const std::vector<CompressedSubgraphs>& subgraphs) {
    uint64_t count = 0;
    for (auto& set : subgraphs) {
      count += set.getNumIsomorphicSubgraphs();
    }
    return count;
  }

  uint64_t getNumSubgraphs(const std::vector<CompressedSubgraphs>& subgraphs) {
    uint64_t count = 0;
    for (auto& set : subgraphs) {
      count += set.getNumSubgraphs();
    }
    return count;
  }

  void batchDFSExecute(const Graph* g, ExecutionPlan* plan) {
    LOG(INFO) << FLAGS_num_cores << " threads";
    ThreadPool threads(FLAGS_num_cores, plan);
    auto& seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    if (plan->isInCover(plan->getRootQueryVertexID())) {
      for (size_t i = 0; i < seeds.size(); i += BATCH_SIZE) {
        size_t end = std::min(i + BATCH_SIZE, seeds.size());
        threads.addInitTask(0, std::vector<CompressedSubgraphs>(seeds.begin() + i, seeds.begin() + end), g);
      }
    } else {
      threads.addInitTask(
          0, std::vector<CompressedSubgraphs>({CompressedSubgraphs(std::make_shared<std::vector<VertexID>>(seeds))}),
          g);
    }
    threads.start();
  }

  void batchDFSExecuteST(const Graph* g, ExecutionPlan* plan) {
    auto seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    if (plan->isInCover(plan->getRootQueryVertexID())) {
      plan->getOperators().handleInput(g, std::vector<CompressedSubgraphs>(seeds.begin(), seeds.end()));
    } else {
      std::vector<CompressedSubgraphs> input;
      input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(seeds)));
      plan->getOperators().handleInput(g, input);
    }
  }

  void bfsExecute(const Graph* g, const ExecutionPlan* plan) {
    auto seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    std::vector<CompressedSubgraphs> outputs;
    if (plan->isInCover(plan->getRootQueryVertexID())) {
      std::vector<CompressedSubgraphs> input(seeds.begin(), seeds.end());
      auto current_op = plan->getOperators().root();
      auto op = dynamic_cast<circinus::TraverseOperator*>(current_op);
      uint32_t op_idx = 0;
      while (op != nullptr) {
        op->input(input, g);
        while (op->expand(&outputs, FLAGS_batch_size) > 0) {
        }
        input.clear();
        LOG(INFO) << op_idx++ << ": # groups " << outputs.size() << " # matches " << getNumIsomorphicSubgraphs(outputs)
                  << '/' << getNumSubgraphs(outputs) << " " << op->toString();
        input.swap(outputs);
        current_op = current_op->getNext();
        op = dynamic_cast<circinus::TraverseOperator*>(current_op);
      }
      auto output_op = dynamic_cast<circinus::OutputOperator*>(current_op);
      CHECK(output_op != nullptr);
      output_op->validateAndOutput(input, 0);
    } else {
      LOG(ERROR) << "this branch should not be reached";
    }
  }

  std::vector<QueryVertexID> getOrder(const std::string& order_str, uint32_t size) {
    std::vector<QueryVertexID> order;
    if (order_str.empty()) return order;
    order.reserve(size);
    QueryVertexID v = 0;
    for (auto c : order_str) {
      if (c == ' ') {
        order.push_back(v);
        v = 0;
      } else {
        v = v * 10 + (c - '0');
      }
    }
    order.push_back(v);
    CHECK_EQ(order.size(), size);
    return order;
  }

  template <VertexCoverStrategy vcs>
  void run(const std::string& graph_path, const std::string& query_path, std::ostream* out) {
    auto start_loading = std::chrono::steady_clock::now();
    Graph g(FLAGS_data_dir + "/" + graph_path);  // load data graph
    auto end_loading = std::chrono::steady_clock::now();
    LOG(INFO) << "========================";
    LOG(INFO) << "graph " << graph_path << " query " << query_path;
    QueryGraph q(FLAGS_data_dir + "/" + query_path);  // load query graph
    auto use_order = getOrder(FLAGS_match_order, q.getNumVertices());
    auto start_filter = std::chrono::steady_clock::now();
    auto candidates = getCandidateSets(g, q);  // get candidates for each query vertex
    auto end_filter = std::chrono::steady_clock::now();
    std::vector<double> candidate_cardinality;
    candidate_cardinality.reserve(candidates.size());
    for (auto& set : candidates) {
      double size = set.size();
      candidate_cardinality.push_back(std::log2(size));
    }
    Profiler profiler;
    NaivePlanner planner(&q, &candidate_cardinality);
    ExecutionPlan* plan;
    if (vcs == Static) {
      plan = planner.generatePlan(use_order);
    } else if (vcs == Eager) {
      plan = planner.generatePlanWithEagerDynamicCover(use_order);
    } else {
      LOG(ERROR) << "Unknown vertex cover strategy " << FLAGS_vertex_cover;
      return;
    }
    plan->setCandidateSets(candidates);  // swap
    plan->printPhysicalPlan();
    // plan->printLabelFrequency();
    plan->getOutputs().init(FLAGS_num_cores).limit(FLAGS_match_limit);
    LOG(INFO) << "limit per thread " << plan->getOutputs().getLimitPerThread();
    LOG(INFO) << &profiler;
    auto start_execution = std::chrono::steady_clock::now();
    // ProfilerStart("benchmark.prof");
    // bfsExecute(&g, plan);
    if (FLAGS_num_cores == 1) {
      LOG(INFO) << "batchDFSExecuteST";
      batchDFSExecuteST(&g, plan);
    } else {
      batchDFSExecute(&g, plan);
    }
    // ProfilerStop();
    auto n_matches = plan->getOutputs().getCount();
    auto end = std::chrono::steady_clock::now();
    auto& order = planner.getMatchingOrder();
    std::stringstream ss;
    for (auto v : order) {
      ss << v << ' ';
    }

    if (FLAGS_profile_file != "") {
      std::ofstream profile_stream;
      profile_stream.open(FLAGS_profile_file);
      profiler.profile(&profile_stream);
      profile_stream.close();
    }

    (*out) << toSeconds(start_loading, end_loading) << ',' << toSeconds(start_filter, end_filter) << ','
           << toSeconds(end_filter, start_execution) << ',' << toSeconds(start_execution, end) << ',' << n_matches
           << ',' << ss.str() << '\n';
  }
};

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_log_dir == "") {
    google::LogToStderr();
  }

  Benchmark benchmark;
  // header: dataset,query_size,query_mode,query_index,load_time,filter_time,plan_time,enumerate_time,n_embeddings,order
  std::ostream* out;
  std::ofstream fstream;
  if (FLAGS_output_file != "") {
    fstream.open(FLAGS_output_file);
    CHECK(fstream.is_open());
    out = &fstream;
  } else {
    out = &std::cout;
  }
  if (FLAGS_vertex_cover == "static") {
    benchmark.run<Static>(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  } else {
    benchmark.run<Eager>(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  }
  fstream.close();

  return 0;
}

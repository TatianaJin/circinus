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
#include <random>
#include <sstream>
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
#include "ops/types.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/hashmap.h"
#include "utils/profiler.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::GraphType;
using circinus::LDFScan;
using circinus::NaivePlanner;
using circinus::NLFFilter;
using circinus::CFLFilter;
using circinus::CFLOrder;
using circinus::DPISOFilter;
using circinus::OrderBase;
using circinus::TSOOrder;
using circinus::TSOFilter;
using circinus::GQLFilter;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::Task;
using circinus::ThreadPool;
using circinus::VertexID;
using circinus::Profiler;
using circinus::DPISOFilter;
using circinus::TSOOrder;
using circinus::TSOFilter;
using circinus::GQLFilter;
using circinus::OrderBase;
using circinus::CoverNode;
using circinus::QueryType;
using circinus::TraverseOperator;

#define BATCH_SIZE FLAGS_batch_size
#define toSeconds(start, end) \
  (((double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1e9)

DEFINE_string(output_file, "", "The output file path");
DEFINE_string(dataset, "dblp", "The dataset to use");
DEFINE_string(query_mode, "dense", "Dense or sparse query");
DEFINE_uint64(query_size, 8, "The query size.");
DEFINE_uint64(match_limit, 1e5, "The limit of matches to find");
DEFINE_uint64(query_index, 1, "The index of query in the same category");
DEFINE_string(match_order, "", "Matching order");
DEFINE_string(filter, "nlf", "filter");
DEFINE_string(profile_file, "", "profile file");
DEFINE_string(vertex_cover, "static", "Vertex cover strategy: static, eager, all");
DEFINE_string(batch_file, "", "Batch query file");
DEFINE_bool(bipartite_graph, false, "use bipartite graph or not");

enum VertexCoverStrategy : uint32_t { Static = 0, Eager, Sample, Dynamic, All };

class Benchmark {
 protected:
  const std::vector<std::string> datasets_ = {"dblp",    "eu2005",  "hprd",  "human",
                                              "patents", "wordnet", "yeast", "youtube"};

 public:
  template <VertexCoverStrategy vcs>
  void run(const std::string& dataset, uint32_t query_size, const std::string& query_mode, uint32_t index,
           std::ostream* out) {
    auto query_path = dataset + "/query_graph/query_" + query_mode + "_" + std::to_string(query_size) + "_" +
                      std::to_string(index) + ".graph";

    (*out) << dataset << ',' << query_size << ',' << query_mode << ',' << index << ',';
    run<vcs>(dataset, query_path, out);
  }

  template <VertexCoverStrategy vcs>
  void run(double load_time, const Graph& data_graph, const std::string& dataset, uint32_t query_size,
           const std::string& query_mode, uint32_t index, std::ostream* out) {
    auto query_path = dataset + "/query_graph/query_" + query_mode + "_" + std::to_string(query_size) + "_" +
                      std::to_string(index) + ".graph";

    (*out) << dataset << ',' << query_size << ',' << query_mode << ',' << index << ',';
    run<vcs>(load_time, data_graph, query_path, out);
  }

  void loadDataset(Graph& data_graph, const std::string& dataset, double& load_time) {
    auto graph_path = dataset + "/data_graph/" + dataset + ".graph";
    std::ifstream input(FLAGS_data_dir + "/" + graph_path + ".bin", std::ios::binary);
    if (input.is_open()) {
      auto start_loading = std::chrono::steady_clock::now();
      data_graph.loadCompressed(input);
      auto end_loading = std::chrono::steady_clock::now();
      LOG(INFO) << "========================";
      load_time = toSeconds(start_loading, end_loading);
      LOG(INFO) << "graph " << graph_path << ".bin load time " << load_time;
      return;
    }
    auto start_loading = std::chrono::steady_clock::now();
    data_graph = Graph(FLAGS_data_dir + "/" + graph_path);
    auto end_loading = std::chrono::steady_clock::now();
    data_graph.saveAsBinary(FLAGS_data_dir + "/" + graph_path + ".bin");
    LOG(INFO) << "========================";
    load_time = toSeconds(start_loading, end_loading);
    LOG(INFO) << "graph " << graph_path << " load time " << load_time;
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
        if (FLAGS_filter == "dpiso" || FLAGS_filter == "ldf") {
          candidates[v].insert(candidates[v].end(), buffer.begin(), buffer.end());
        } else {
          filter.Filter(g, buffer, &candidates[v]);
        }
        buffer.clear();
      }
      candidate_size[v] = candidates[v].size();
    }
    if (FLAGS_filter == "cfl") {
      CFLOrder cfl_order;
      QueryVertexID start_vertex = cfl_order.getStartVertex(&g, &q, candidate_size);
      // LOG(INFO) << "cfl order get start vertex " << start_vertex;
      CFLFilter cfl_filter(&q, &g, start_vertex);
      cfl_filter.Filter(candidates);
    } else if (FLAGS_filter == "dpiso") {
      OrderBase dpiso_order;
      QueryVertexID start_vertex = dpiso_order.getStartVertex(&g, &q, candidate_size);
      // LOG(INFO) << "dpiso order get start vertex " << start_vertex;
      DPISOFilter dpiso_filter(&q, &g, start_vertex);
      dpiso_filter.Filter(candidates);
    } else if (FLAGS_filter == "tso") {
      TSOOrder tso_order;
      QueryVertexID start_vertex = tso_order.getStartVertex(&g, &q, candidate_size);
      // LOG(INFO) << "tso order get start vertex " << start_vertex;
      TSOFilter tso_filter(&q, &g, start_vertex);
      tso_filter.Filter(candidates);
    } else if (FLAGS_filter == "gql") {
      std::vector<std::vector<VertexID>> gql_candidates;
      circinus::unordered_map<QueryVertexID, circinus::unordered_set<VertexID>> gql_map;
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        circinus::unordered_set<VertexID> gql_set;
        for (auto i : candidates[v]) {
          gql_set.insert(i);
        }
        gql_map[v] = std::move(gql_set);
      }
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        GQLFilter gql_filter(&q, v, &gql_map);
        gql_filter.preFilter(g, candidates[v]);
      }
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        GQLFilter gql_filter(&q, v, &gql_map);
        std::vector<VertexID> tempvec;
        gql_filter.Filter(g, candidates[v], &tempvec);
        gql_candidates.push_back(std::move(tempvec));
      }
      return gql_candidates;
    }

    //  for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
    //   LOG(INFO) << "vertex " << v << " " << candidate_size[v] << "/" << candidates[v].size();
    //}
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

  void batchDFSExecute(const Graph* g, ExecutionPlan* plan, Profiler* profiler) {
    LOG(INFO) << FLAGS_num_cores << " threads";
    ThreadPool threads(FLAGS_num_cores, plan);
    auto& seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    plan->setProfiler(profiler);
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

  std::vector<circinus::GraphView<circinus::BipartiteGraph>*> setupBipartiteGraphs(const Graph* g,
                                                                                   ExecutionPlan* plan) {
    std::vector<circinus::GraphView<circinus::BipartiteGraph>*> data_graphs_for_operators;
    auto& opTree = plan->getOperators();
    auto candidate_sets = plan->getCandidateSets();
    size_t len = opTree.getOperatorSize() - 1;
    data_graphs_for_operators.reserve(len);
    for (size_t i = 0; i < len; ++i) {
      auto op = opTree.getOperator(i);
      auto traverse_op = dynamic_cast<TraverseOperator*>(op);
      data_graphs_for_operators.emplace_back(
          new circinus::GraphView<circinus::BipartiteGraph>(traverse_op->computeBipartiteGraphs(g, *candidate_sets)));
    }
    return data_graphs_for_operators;
  }

  void batchDFSExecuteST(const Graph* g, ExecutionPlan* plan) {
    auto seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    if (FLAGS_bipartite_graph) {
      auto data_graphs_for_operators = setupBipartiteGraphs(g, plan);
      if (plan->isInCover(plan->getRootQueryVertexID()) &&
          (FLAGS_vertex_cover == "static" || FLAGS_vertex_cover == "all" ||
           plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0)) {
        plan->getOperators().execute(data_graphs_for_operators,
                                     std::vector<CompressedSubgraphs>(seeds.begin(), seeds.end()), 0);
      } else {
        std::vector<CompressedSubgraphs> input;
        input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(seeds)));
        plan->getOperators().execute(data_graphs_for_operators, input, 0);
      }
      for (auto& dg : data_graphs_for_operators) {
        delete dg;
      }
      return;
    }
    if (plan->isInCover(plan->getRootQueryVertexID()) &&
        (FLAGS_vertex_cover == "static" || FLAGS_vertex_cover == "all" ||
         plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0)) {
      plan->getOperators().execute(g, std::vector<CompressedSubgraphs>(seeds.begin(), seeds.end()), 0);
    } else {
      std::vector<CompressedSubgraphs> input;
      input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(seeds)));
      plan->getOperators().execute(g, input, 0);
    }
  }

  void batchDFSProfileST(const Graph* g, ExecutionPlan* plan) {
    auto seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    if (FLAGS_bipartite_graph) {
      auto data_graphs_for_operators = setupBipartiteGraphs(g, plan);
      if (plan->isInCover(plan->getRootQueryVertexID()) &&
          (FLAGS_vertex_cover == "static" || FLAGS_vertex_cover == "all" ||
           plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0)) {
        plan->getOperators().profile(data_graphs_for_operators,
                                     std::vector<CompressedSubgraphs>(seeds.begin(), seeds.end()), 0);
      } else {
        std::vector<CompressedSubgraphs> input;
        input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(seeds)));
        plan->getOperators().profile(data_graphs_for_operators, input, 0);
      }
      return;
    }
    if (plan->isInCover(plan->getRootQueryVertexID()) &&
        (FLAGS_vertex_cover == "static" || FLAGS_vertex_cover == "all" ||
         plan->getToKeyLevel(plan->getRootQueryVertexID()) == 0)) {
      plan->getOperators().profile(g, std::vector<CompressedSubgraphs>(seeds.begin(), seeds.end()), FLAGS_profile, 0);
    } else {
      std::vector<CompressedSubgraphs> input;
      input.emplace_back(std::make_shared<std::vector<VertexID>>(std::move(seeds)));
      plan->getOperators().profile(g, input, FLAGS_profile, 0);
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
        auto start = std::chrono::steady_clock::now();
        op->input(input, g);
        while (op->expand(&outputs, FLAGS_batch_size) > 0) {
        }
        auto end = std::chrono::steady_clock::now();
        input.clear();
        LOG(INFO) << op_idx++ << ": # groups " << outputs.size() << " # matches " << getNumIsomorphicSubgraphs(outputs)
                  << '/' << getNumSubgraphs(outputs) << " " << op->toString() << " " << toSeconds(start, end) << "s";
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

  void sampleBfsExecute(const Graph* g, const ExecutionPlan* plan, const std::vector<uint32_t>& order_idx,
                        std::vector<std::vector<double>>& cardinality, std::vector<double>& level_cost) {
    auto seeds = plan->getCandidateSet(plan->getRootQueryVertexID());
    std::vector<CompressedSubgraphs> outputs;
    cardinality.resize(order_idx.size());
    std::vector<double> sample_ratio;

    const uint32_t sample_num = 5000;
    std::vector<VertexID> sample_seeds;
    std::random_device rd;
    std::sample(seeds.begin(), seeds.end(), std::back_inserter(sample_seeds), sample_num, std::mt19937(rd()));
    std::vector<CompressedSubgraphs> input(sample_seeds.begin(), sample_seeds.end());
    auto current_op = plan->getOperators().root();
    auto op = dynamic_cast<circinus::TraverseOperator*>(current_op);
    uint32_t op_idx = 0;
    sample_ratio.emplace_back(seeds.size() / (double)sample_seeds.size());
    // the first level
    cardinality[0].emplace_back(seeds.size());
    while (op != nullptr) {
      auto start = std::chrono::steady_clock::now();
      op->input(input, g);
      while (op->expandAndProfile(&outputs, FLAGS_batch_size, 1) > 0) {
      }
      auto end = std::chrono::steady_clock::now();
      double time = toSeconds(start, end);
      level_cost.emplace_back(time);
      LOG(INFO) << time << "  " << op->getIntersectionCount() << "  " << op->getTotalIntersectionInputSize() << "  "
                << op->getTotalIntersectionOutputSize() << "  " << op->getDistinctIntersectionCount();

      LOG(INFO) << op_idx++ << ": # groups " << outputs.size() << " # matches " << getNumIsomorphicSubgraphs(outputs)
                << '/' << getNumSubgraphs(outputs) << " " << op->toString() << " " << toSeconds(start, end) << "s";
      cardinality[op_idx].resize(order_idx.size());
      for (uint32_t i = 0; i <= op_idx; ++i) {
        circinus::unordered_set<VertexID> vertex_set;
        for (const auto& output : outputs) {
          vertex_set.insert(output.getKeyVal(i));
        }
        cardinality[op_idx][i] = vertex_set.size();
      }
      input.clear();
      std::sample(outputs.begin(), outputs.end(), std::back_inserter(input), sample_num, std::mt19937(rd()));

      sample_ratio.emplace_back(outputs.size() / (double)input.size());
      outputs.clear();
      current_op = current_op->getNext();
      op = dynamic_cast<circinus::TraverseOperator*>(current_op);
    }

    for (uint32_t i = 0; i < cardinality.size(); ++i) {
      std::string s = "level " + std::to_string(i) + " ";
      for (uint32_t j = 0; j < i + 1; ++j) {
        s += std::to_string(cardinality[i][j]) + " ";
      }
      LOG(INFO) << s;
    }
    for (uint32_t i = 0; i < level_cost.size(); ++i) {
      DLOG(INFO) << "op " << i << " time usage " << level_cost[i];
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

  void run_sample(const Graph& g, QueryGraph& q, std::vector<std::vector<VertexID>>& candidates,
                  const std::vector<QueryVertexID>& use_order, std::vector<std::vector<double>>& cardinality,
                  std::vector<double>& level_cost) {
    std::vector<double> candidate_cardinality;
    candidate_cardinality.reserve(candidates.size());
    for (auto& set : candidates) {
      double size = set.size();
      candidate_cardinality.push_back(size);
    }

    NaivePlanner planner(&q, &candidate_cardinality);
    ExecutionPlan* plan;
    plan = planner.generatePlanWithoutCompression(use_order);
    plan->setCandidateSets(candidates);  // swap
    plan->printPhysicalPlan();
    plan->getOutputs().init(FLAGS_num_cores).limit(FLAGS_match_limit);

    // order to order_idx
    std::vector<uint32_t> order_idx(use_order.size(), 0);
    for (uint32_t i = 0; i < use_order.size(); ++i) {
      order_idx[use_order[i]] = i;
    }
    auto start = std::chrono::steady_clock::now();
    std::vector<std::vector<double>> cardinality_tmp;
    sampleBfsExecute(&g, plan, order_idx, cardinality_tmp, level_cost);
    cardinality.resize(cardinality_tmp.size());
    for (uint32_t i = 0; i < cardinality_tmp.size(); ++i) {
      cardinality[i].resize(use_order.size());
      for (uint32_t j = 0; j < cardinality_tmp[i].size(); ++j) {
        cardinality[i][use_order[j]] = cardinality_tmp[i][j];
      }
    }
    auto end = std::chrono::steady_clock::now();
    LOG(INFO) << "sample execution time " << toSeconds(start, end);
  }

  template <VertexCoverStrategy vcs>
  void run(const std::string& dataset, const std::string& query_path, std::ostream* out) {
    Graph g;
    double load_time;
    loadDataset(g, dataset, load_time);
    LOG(INFO) << "query " << query_path;
    run<vcs>(load_time, g, query_path, out);
  }

  template <VertexCoverStrategy vcs>
  void run(double load_time, const Graph& g, const std::string& query_path, std::ostream* out) {
    QueryGraph q(FLAGS_data_dir + "/" + query_path);  // load query graph
    auto use_order = getOrder(FLAGS_match_order, q.getNumVertices());
    auto start_filter = std::chrono::steady_clock::now();
    auto candidates = getCandidateSets(g, q);  // get candidates for each query vertex
    auto end_filter = std::chrono::steady_clock::now();
    std::vector<double> candidate_cardinality;
    candidate_cardinality.reserve(candidates.size());
    for (auto& set : candidates) {
      double size = set.size();
      candidate_cardinality.push_back(size);
    }
    Profiler profiler;  // TODO(tatiana): deprecate
    NaivePlanner planner(&q, &candidate_cardinality,
                         FLAGS_bipartite_graph ? GraphType::BipartiteGraphView : GraphType::Normal);
    ExecutionPlan* plan;
    if (vcs == Static) {
      plan = planner.generatePlan(use_order);
    } else if (vcs == Eager) {
      plan = planner.generatePlanWithEagerDynamicCover(use_order);
    } else if (vcs == All) {
      plan = planner.generatePlanWithoutCompression(use_order);
    } else if (vcs == Sample) {
      planner.generateOrder(use_order);
      const auto& match_order = planner.getMatchingOrder();
      std::vector<std::vector<double>> cardinality;
      std::vector<double> level_cost;
      run_sample(g, q, candidates, match_order, cardinality, level_cost);
      planner.setCandidateCardinality(&cardinality.back());
      planner.generateCoverNode(cardinality);
      plan = planner.generatePlanWithSampleExecution(cardinality, level_cost);
    } else if (vcs == Dynamic) {
      planner.generateOrder(use_order);
      std::vector<std::vector<double>> cardinality(candidate_cardinality.size(), candidate_cardinality);
      planner.generateCoverNode(cardinality);
      plan = planner.generatePlanWithDynamicCover();
    } else {
      LOG(ERROR) << "Unknown vertex cover strategy " << FLAGS_vertex_cover;
      return;
    }
    LOG(INFO) << "Generate Plan Finish.";
    plan->setCandidateSets(candidates);  // swap
    plan->printPhysicalPlan();
    // plan->printLabelFrequency();
    plan->getOutputs().init(FLAGS_num_cores).limit(FLAGS_match_limit);
    LOG(INFO) << "limit per thread " << plan->getOutputs().getLimitPerThread();
    auto start_execution = std::chrono::steady_clock::now();
    // ProfilerStart("benchmark.prof");
    if (FLAGS_num_cores == 1) {
      if (FLAGS_profile > 0) {
        batchDFSProfileST(&g, plan);
      } else {
        LOG(INFO) << "batchDFSExecuteST";
        batchDFSExecuteST(&g, plan);
      }
      // bfsExecute(&g, plan);
    } else {
      batchDFSExecute(&g, plan, &profiler);
    }
    // ProfilerStop();
    auto n_matches = plan->getOutputs().getCount();
    auto end = std::chrono::steady_clock::now();
    auto& order = planner.getMatchingOrder();
    std::stringstream ss;
    for (auto v : order) {
      ss << v << ' ';
    }

    if (FLAGS_profile_file != "" && FLAGS_num_cores > 1) {
      std::ofstream profile_stream;
      profile_stream.open(FLAGS_profile_file);
      profiler.profile(&profile_stream);
      profile_stream.close();
    } else if (FLAGS_profile_file != "") {
      std::ofstream profile_stream;
      profile_stream.open(FLAGS_profile_file);
      plan->printProfiledPlan(profile_stream);
      profile_stream.close();
    }

    (*out) << load_time << ',' << toSeconds(start_filter, end_filter) << ',' << toSeconds(end_filter, start_execution)
           << ',' << toSeconds(start_execution, end) << ',' << n_matches << ',' << ss.str();
  }
};

class QueryConfig {
 public:
  std::string dataset;
  uint32_t query_size;
  std::string query_mode;
  uint32_t query_index;
  std::string match_order;

  bool skipConfig() const { return skip_; }

  void readNextConfig(std::istream& str) {
    skip_ = false;
    std::string::size_type pos = 0;
    auto old_pos = pos;
    std::getline(str, line_);
    // parse dataset
    if ((pos = line_.find(',', pos)) == std::string::npos) {
      skip_ = true;
      return;
    }
    dataset = std::string(&line_[old_pos], pos - old_pos);
    if (dataset.substr(0, 7) == "dataset") {
      skip_ = true;
      return;
    }
    old_pos = ++pos;
    // parse query_size
    CHECK_NE((pos = line_.find(',', pos)), std::string::npos);
    line_[pos] = '\0';
    query_size = std::atoi(&line_[old_pos]);
    old_pos = ++pos;
    // parse query_mode
    CHECK_NE((pos = line_.find(',', pos)), std::string::npos);
    query_mode = std::string(&line_[old_pos], pos - old_pos);
    old_pos = ++pos;
    // parse query_index
    CHECK_NE((pos = line_.find(',', pos)), std::string::npos);
    line_[pos] = '\0';
    query_index = std::atoi(&line_[old_pos]);
    old_pos = ++pos;
    // parse match_order
    pos = line_.find(',', pos);
    if (pos == line_.npos) {
      pos = line_.size();
    }
    if (line_[pos - 1] == ' ') {
      match_order = std::string(&line_[old_pos], pos - old_pos - 1);
    } else {
      match_order = std::string(&line_[old_pos], pos - old_pos);
    }
  }

  std::ostream& toString(std::ostream& ss) {
    if (skip_) return ss;
    ss << "-dataset " << dataset << " -query_size " << query_size << " -query_mode " << query_mode << " -query_index "
       << query_index << " -match_order " << match_order << std::endl;
    return ss;
  }

 private:
  std::string line_;
  bool skip_ = false;
};

std::istream& operator>>(std::istream& str, QueryConfig& config) {
  config.readNextConfig(str);
  return str;
}

void run_benchmark(const std::string& query_file, std::ostream* out) {
  std::ifstream query_f(query_file);
  LOG(INFO) << "batch query file" << query_file;
  CHECK(query_f.is_open());

  (*out) << "dataset,query_size,query_mode,query_index,load_time,filter_time,plan_time,enumerate_time,n_embeddings,"
            "order\n";
  Benchmark benchmark;
  QueryConfig config;
  std::string dataset = "";
  Graph data_graph;
  double load_time = 0;
  std::string profile_prefix = "/data/share/users/byli/circinus/evaluation/profile/";
  while (query_f >> config) {
    if (config.skipConfig()) continue;
    config.toString(LOG(INFO) << ">>>>>>>>>>>>>>>>> query config -vertex_cover " << FLAGS_vertex_cover << " -filter "
                              << FLAGS_filter << " -match_limit " << FLAGS_match_limit << ' ');
    // load graph if not cached
    FLAGS_profile_file = profile_prefix + config.dataset + '_' + config.query_mode + '_' +
                         std::to_string(config.query_size) + '_' + std::to_string(config.query_index) + '_' +
                         FLAGS_filter;
    LOG(INFO) << "-------------" << FLAGS_profile_file;
    if (config.dataset != dataset) {
      dataset = config.dataset;
      benchmark.loadDataset(data_graph, dataset, load_time);
    }
    FLAGS_match_order = config.match_order;
    if (FLAGS_vertex_cover == "static") {
      benchmark.run<Static>(load_time, data_graph, config.dataset, config.query_size, config.query_mode,
                            config.query_index, out);
    } else if (FLAGS_vertex_cover == "all") {
      benchmark.run<All>(load_time, data_graph, config.dataset, config.query_size, config.query_mode,
                         config.query_index, out);
    } else if (FLAGS_vertex_cover == "eager") {
      benchmark.run<Eager>(load_time, data_graph, config.dataset, config.query_size, config.query_mode,
                           config.query_index, out);
    } else if (FLAGS_vertex_cover == "dynamic") {
      benchmark.run<Dynamic>(load_time, data_graph, config.dataset, config.query_size, config.query_mode,
                             config.query_index, out);
    } else {
      benchmark.run<Sample>(load_time, data_graph, config.dataset, config.query_size, config.query_mode,
                            config.query_index, out);
    }
    (*out) << std::endl;
  }
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

  Benchmark benchmark;
  // header: dataset,query_size,query_mode,query_index,load_time,filter_time,plan_time,enumerate_time,n_embeddings,order
  std::ostream* out;
  std::ofstream fstream;
  if (FLAGS_output_file != "") {
    fstream.open(FLAGS_output_file, std::ios::app);
    CHECK(fstream.is_open());
    out = &fstream;
  } else {
    out = &std::cout;
  }
  if (FLAGS_match_limit == 0) {
    FLAGS_match_limit = ~0ull;
  }
  // FLAGS_profile = (FLAGS_profile_file != "");

  if (FLAGS_batch_file != "") {
    run_benchmark(FLAGS_batch_file, out);
    return 0;
  }

  if (FLAGS_vertex_cover == "static") {
    benchmark.run<Static>(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  } else if (FLAGS_vertex_cover == "all") {
    benchmark.run<All>(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  } else if (FLAGS_vertex_cover == "eager") {
    benchmark.run<Eager>(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  } else if (FLAGS_vertex_cover == "dynamic") {
    benchmark.run<Dynamic>(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  } else {
    benchmark.run<Sample>(FLAGS_dataset, FLAGS_query_size, FLAGS_query_mode, FLAGS_query_index, out);
  }
  (*out) << std::endl;

  fstream.close();

  return 0;
}

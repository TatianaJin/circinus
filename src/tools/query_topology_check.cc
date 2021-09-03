// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
// the License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

#include <string>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "algorithms/k_clique.h"
#include "algorithms/k_cycle.h"
#include "algorithms/k_star.h"
#include "graph/query_graph.h"
#include "utils/file_utils.h"

using circinus::KClique;
using circinus::KCycle;
using circinus::KStar;
using circinus::Path;
using circinus::QueryGraph;

template <uint32_t K>
inline void run(int argc, char** argv) {
  // log header
  std::stringstream header;
  header << "query";
  for (uint32_t k = 3; k <= K; ++k) header << ',' << k << "clique";
  for (uint32_t k = 4; k <= K; ++k) header << ',' << k << "cycle";
  for (uint32_t k = 2; k <= K; ++k) header << ',' << k << "star";
  header << std::endl;
  std::cout << header.str();
  // for each query graph, compute the stats
  for (int i = 2; i < argc; ++i) {
    QueryGraph q(argv[i]);
    KClique<K> clique_solver(&q);
    KCycle<K> cycle_solver(&q);
    KStar<K> star_solver(&q);
    std::stringstream ss;
    ss << argv[i];
    CHECK_EQ(clique_solver.getCliqueCount(3), cycle_solver.getCycleCount(3));
    for (uint32_t k = 3; k <= K; ++k) {
      ss << ',' << clique_solver.getCliqueCount(k);
    }
    for (uint32_t k = 4; k <= K; ++k) {
      ss << ',' << cycle_solver.getCycleCount(k);
    }
    for (uint32_t k = 2; k <= K; ++k) {
      ss << ',' << star_solver.getStarCount(k);
    }
    ss << std::endl;
    std::cout << ss.str();
  }
}

template <uint32_t K>
inline void label_count(int argc, char** argv) {
  // for each query graph, compute the stats
  for (int i = 2; i < argc; ++i) {
    QueryGraph q(argv[i]);
    KClique<K, true> clique_solver(&q);
    KCycle<K, true> cycle_solver(&q);
    std::stringstream ss;
    ss << argv[i];
    for (auto& pair : clique_solver.getLabeledCliqueCounts()) {
      ss << ',' << pair.first << ':' << pair.second;
    }
    for (auto& pair : cycle_solver.getLabeledCycleCounts()) {
      ss << ',' << pair.first << ':' << pair.second;
    }
    for (circinus::QueryVertexID v = 0; v < q.getNumVertices(); ++v) {
      auto nbrs = q.getOutNeighbors(v);
      ss << ",nbr";
      std::vector<circinus::LabelID> labels(nbrs.second);
      for (uint32_t nb = 0; nb < nbrs.second; ++nb) {
        labels[nb] = q.getVertexLabel(nbrs.first[nb]);
      }
      std::sort(labels.begin(), labels.end());
      for (auto label : labels) ss << '|' << label;
      ss << ":1";
    }
    ss << std::endl;
    std::cout << ss.str();
  }
}

template <uint32_t i>
void unroll(uint32_t n, int argc, char** argv) {
  if (i == n) return label_count<i>(argc, argv);
  unroll<i + 1>(n, argc, argv);
}

template <>
void unroll<33>(uint32_t n, int argc, char** argv) {
  LOG(WARNING) << "Expecting n <= 32, got n = " << n;
}

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  google::SetUsageMessage("./QueryTopologyCheck QUERY_GRAPH_SIZE QUERY_GRAPH_PATHS ...");
  if (argc < 3) {
    std::cout << gflags::ProgramUsage() << std::endl;
    return 0;
  }
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_log_dir == "") google::LogToStderr();

  uint32_t query_graph_size = std::atoi(argv[1]);
  CHECK_GE(query_graph_size, 4);
  unroll<4>(query_graph_size, argc, argv);
}

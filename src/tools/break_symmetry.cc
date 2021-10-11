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

#include "algorithms/automorphism_check.h"
#include "graph/query_graph.h"

using circinus::AutomorphismCheck;
using circinus::QueryGraph;
using circinus::QueryVertexID;

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_log_dir == "") google::LogToStderr();

  if (argc == 1) {
    std::cout << "usage: ./BreakSymmetry QUERY_GRAPH_PATH" << std::endl;
    return 0;
  }

  std::string query_path = argv[1];

  QueryGraph q(query_path);

  AutomorphismCheck ac(q);
  auto conditions = ac.getPartialOrder();
  std::cout << "Partial orders:" << std::endl;
  for (auto& p : conditions) {
    std::cout << p.first << ' ' << p.second << std::endl;
  }
  return 0;
}

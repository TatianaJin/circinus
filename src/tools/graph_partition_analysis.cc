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

#include "glog/logging.h"

#include "graph/partitioned_graph.h"
#include "utils/file_utils.h"

using circinus::ReorderedPartitionedGraph;

DEFINE_string(graph, "", "The binary file path of the partitioned graph.");

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  google::SetUsageMessage("./GraphPartitionAnalysis --graph GRAPH_PATH");

  if (FLAGS_log_dir == "") google::LogToStderr();

  if (FLAGS_graph.empty()) {
    std::cout << gflags::ProgramUsage() << std::endl;
    return 0;
  }

  auto input = circinus::openFile(FLAGS_graph);
  auto graph = circinus::GraphBase::loadGraphFromBinary(input);
  auto& pg = *dynamic_cast<ReorderedPartitionedGraph*>(graph.get());

  pg.showPartitionInfo();
  return 0;
}

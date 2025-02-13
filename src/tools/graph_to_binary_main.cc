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

#include <iostream>

#include "glog/logging.h"

#include "graph/graph.h"
#include "graph/partitioned_graph.h"

using circinus::Graph;
using circinus::ReorderedPartitionedGraph;

DEFINE_int32(partition, 0,
             "By default partition=0, i.e., no partitioning or reordering. If 1, graph is reordered without "
             "partitioning. Otherwise, graph is partitioned into N partitions");
DEFINE_bool(sort_by_degree, true,
            "By default sort_by_degree=true, that is, when reordering is enabled, degree sorting is applied as the "
            "minor order.");

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  google::SetUsageMessage(
      "./GraphToBinary graph_path [output_path] [--partition=N] [--sort_by_degree=true]"
      "\n\tpartition\n\t\t"
      "By default partition=0, i.e., no partitioning or reordering. If 1, graph is reordered "
      "without partitioning. Otherwise, graph is partitioned into N partitions."
      "\n\tsort_by_degree\n\t\t"
      "By default sort_by_degree=true, that is, when reordering is enabled, degree sorting "
      "is applied as the minor order.");
  if (FLAGS_log_dir == "") {
    google::LogToStderr();
  }
  if (argc == 1) {
    std::cout << gflags::ProgramUsage() << std::endl;
    return 0;
  }
  std::string input = argv[1];
  std::string output;
  if (argc == 3) {
    output = argv[2];
  } else if (FLAGS_partition == 0) {
    output = input + ".bin";
  } else {
    output = input + "_ordered.bin";
  }
  if (FLAGS_partition == 0) {
    LOG(INFO) << "Transform graph " << input << " to binary, with output path " << output;
    Graph g(input);
    g.saveAsBinary(output);
  } else {
    LOG(INFO) << "Transform graph " << input << " as ReorderedPartitionedGraph to binary, with output path " << output;
    ReorderedPartitionedGraph g(input, FLAGS_partition, FLAGS_sort_by_degree);
    LOG(INFO) << "Loaded graph, now saving as a binary file";
    g.saveAsBinary(output);
    LOG(INFO) << "Saved to " << output;
  }
}

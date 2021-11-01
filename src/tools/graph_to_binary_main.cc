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
DEFINE_string(input_format, "", "Default is our input graph format, option: [snap]");
DEFINE_string(partition_file, "", "Default is none, fennel partition file.");

int main(int argc, char** argv) {
#ifndef NDEBUG
  FLAGS_logbuflevel = -1;  // -1 means don't buffer.
#endif
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  google::SetUsageMessage(
      "./GraphToBinary graph_path [output_path] [--partition=N] [--sort_by_degree=true] [--input_format=snap]"
      "\n\tpartition\n\t\t"
      "By default partition=0, i.e., no partitioning or reordering. If 1, graph is reordered "
      "without partitioning. Otherwise, graph is partitioned into N partitions."
      "\n\tsort_by_degree\n\t\t"
      "By default sort_by_degree=true, that is, when reordering is enabled, degree sorting "
      "is applied as the minor order.\n\t\t"
      "By Default input_format is our input graph format, option: [snap]");
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
  } else if (FLAGS_partition == 1 && FLAGS_sort_by_degree) {
    output = input + ".sorted.bin";
  } else {
    output = input.substr(0, input.size() - 9) + ".circinus.bin.p" + std::to_string(FLAGS_partition);
  }

  if (FLAGS_partition == 0) {
    if (FLAGS_input_format == "snap") {
      Graph g;
      g.loadUndirectedGraphEdgeList(input);
      g.saveAsBinary(output);
    } else if (FLAGS_input_format == "tve") {
      Graph g;
      g.loadTVEUndirectedGraphToGrasperGraph(input);
    } else {
      LOG(INFO) << "Transform graph " << input << " to binary, with output path " << output;
      Graph g(input);
      LOG(INFO) << "Loaded graph, now saving as a binary file";
      g.saveAsBinary(output);
    }
  } else {
    if (FLAGS_input_format == "snap") {
      LOG(INFO) << "Transform snap format graph " << input
                << " as ReorderedPartitionedGraph to binary, with output path " << output;
      Graph g;
      g.loadUndirectedGraphEdgeList(input);
      ReorderedPartitionedGraph rg(g, FLAGS_partition, FLAGS_sort_by_degree);
      rg.saveAsBinary(output);
    } else if (FLAGS_partition_file != "") {
      LOG(INFO) << "Transform graph " << input << " as ReorderedPartitionedGraph to binary with partition file "
                << FLAGS_partition_file << ", with output path " << output;
      ReorderedPartitionedGraph g(input, FLAGS_partition_file, FLAGS_partition, FLAGS_sort_by_degree);
      LOG(INFO) << "Loaded graph, now saving as a binary file";
      g.saveAsBinary(output);
      LOG(INFO) << "Saved to " << output;
    } else if (FLAGS_partition == 1 && FLAGS_sort_by_degree) {
      LOG(INFO) << "Transform graph " << input << " to binary, with output path " << output;
      Graph graph(input);
      graph.reorderByDegree();
      graph.saveAsBinary(output);
      LOG(INFO) << "Saved to " << output;
    } else {
      LOG(INFO) << "Transform graph " << input << " as ReorderedPartitionedGraph to binary, with output path "
                << output;
      ReorderedPartitionedGraph g(input, FLAGS_partition, FLAGS_sort_by_degree);
      LOG(INFO) << "Loaded graph, now saving as a binary file";
      g.saveAsBinary(output);
      LOG(INFO) << "Saved to " << output;
    }
  }
}

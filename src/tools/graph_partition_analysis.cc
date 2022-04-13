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

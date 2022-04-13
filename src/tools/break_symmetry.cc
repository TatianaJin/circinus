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
  auto po = ac.getPartialOrder();
  std::cout << "Partial orders:";
  po.printMinimum(std::cout) << std::endl;
  std::cout << "Full constraints:";
  po.printFull(std::cout) << std::endl;
  return 0;
}

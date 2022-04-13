#pragma once

#include <vector>

#include "exec/task.h"
#include "graph/bipartite_graph.h"
#include "graph/graph_view.h"
#include "plan/execution_plan.h"

namespace circinus {

class BipartiteGraphTraverseTask : public TaskBase {
 public:
  static std::vector<circinus::GraphView<circinus::BipartiteGraph>> setupBipartiteGraphs(
      const Graph* g, const std::vector<CandidateSetView>& candidate_sets, const OperatorTree& op_tree) {
    std::vector<circinus::GraphView<circinus::BipartiteGraph>> data_graphs_for_operators;
    size_t len = op_tree.getOperatorSize() - 1;
    data_graphs_for_operators.reserve(len);
    for (size_t i = 0; i < len; ++i) {
      auto op = op_tree.getOperator(i);
      auto traverse_op = dynamic_cast<TraverseOperator*>(op);
      data_graphs_for_operators.emplace_back(
          new circinus::GraphView<circinus::BipartiteGraph>(traverse_op->computeBipartiteGraphs(g, candidate_sets)));
    }
    return data_graphs_for_operators;
  }
};

}  // namespace circinus

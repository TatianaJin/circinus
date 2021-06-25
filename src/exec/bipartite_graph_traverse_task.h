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

#pragma once

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

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

#include <vector>

#include "graph/query_graph.h"

namespace circinus {

class TwoCoreSolver {
 private:
  std::vector<int> core_table_;
  bool generated=0;

 public:
  inline static bool isInCore(const std::vector<int>& core_table, QueryVertexID v) { return core_table[v] > 1; }

  inline static uint32_t getCoreSize(const std::vector<int>& core_table) {
    uint32_t count = 0;
    for (auto val : core_table) {
      count += (val > 1);
    }
    return count;
  }

  const std::vector<int>& get2CoreTable(const QueryGraph* graph);

  std::vector<QueryVertexID> get2CoreVertices(const QueryGraph* graph);

  inline QueryGraph extract2CoreSubgraph(const QueryGraph* graph) {
    return graph->getInducedSubgraph(get2CoreVertices(graph));
  }
};

}  // namespace circinus

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

#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"

namespace circinus {

class OrderBase {
 public:
  virtual uint32_t getStartVertex(const Graph* data_graph, const QueryGraph* query_graph,
                                  std::vector<uint32_t> candidate_size) {
    double min_score = data_graph->getNumVertices();
    QueryVertexID start_vertex = 0;

    for (QueryVertexID v = 0; v < query_graph->getNumVertices(); ++v) {
      uint32_t degree = query_graph->getVertexOutDegree(v);
      double cur_score = candidate_size[v] / (double)degree;
      if (cur_score < min_score) {
        min_score = cur_score;
        start_vertex = v;
      }
    }

    return start_vertex;
  }
};

}  // namespace circinus

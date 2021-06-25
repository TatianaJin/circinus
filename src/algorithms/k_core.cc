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

#include "algorithms/k_core.h"

#include <vector>
#include "graph/query_graph.h"

namespace circinus {

std::vector<QueryVertexID> TwoCoreSolver::get2CoreVertices(const QueryGraph* graph) {
  std::vector<QueryVertexID> ret;
  auto vertices_count = graph->getNumVertices();
  if (!generated) {
    get2CoreTable(graph);
  }
  for (QueryVertexID i = 0; i < vertices_count; ++i) {
    if (isInCore(core_table_, i)) {
      ret.push_back(i);
    }
  }
  return ret;
}

const std::vector<int>& TwoCoreSolver::get2CoreTable(const QueryGraph* graph) {
  if(generated)return core_table_;
  generated=1;
  auto vertices_count = graph->getNumVertices();
  auto max_degree = graph->getGraphMaxDegree();

  std::vector<QueryVertexID> vertices(vertices_count);  // IDs of vertices sorted by degree.
  std::vector<uint32_t> position(vertices_count);       // The position of vertex i in the sorted array.
  std::vector<uint32_t> degree_bin(max_degree + 1, 0);  // i: the number of vertices with degree i
  std::vector<uint32_t> offsets(max_degree + 1, 0);     // i: the offset of the first vertex with degree i.
  core_table_.resize(vertices_count);

  // populate degree_bin
  for (uint32_t i = 0; i < vertices_count; ++i) {
    int degree = graph->getVertexOutDegree(i);
    core_table_[i] = degree;
    degree_bin[degree] += 1;
  }

  // populate offsets according to degree_bin
  for (uint32_t i = 1; i < max_degree + 1; ++i) {
    offsets[i] = offsets[i - 1] + degree_bin[i - 1];
  }

  // count sort to populate vertices
  // here offsets[i] is used as the position of the next vertex of degreed i
  for (uint32_t i = 0; i < vertices_count; ++i) {
    int degree = graph->getVertexOutDegree(i);
    position[i] = offsets[degree];
    vertices[position[i]] = i;
    offsets[degree] += 1;
  }

  // restore offsets
  for (int i = max_degree; i > 0; --i) {
    offsets[i] = offsets[i - 1];
  }
  offsets.front() = 0;

  for (uint32_t i = 0; i < vertices_count; ++i) {
    int v = vertices[i];

    auto neighbors = graph->getOutNeighbors(v);

    for (uint32_t j = 0; j < neighbors.second; ++j) {
      auto u = neighbors.first[j];

      if (core_table_[u] > core_table_[v]) {
        // Get the position and vertex which is with the same degree
        // and at the start position of vertices array.
        int cur_degree_u = core_table_[u];
        int position_u = position[u];
        int position_w = offsets[cur_degree_u];
        auto w = vertices[position_w];

        if (u != w) {
          // Swap u and w.
          position[u] = position_w;
          position[w] = position_u;
          vertices[position_u] = w;
          vertices[position_w] = u;
        }

        offsets[cur_degree_u] += 1;
        core_table_[u] -= 1;
      }
    }
  }

  return core_table_;
}

}  // namespace circinus

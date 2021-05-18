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

#include <algorithm>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "utils/hashmap.h"

namespace circinus {

struct FilterContext {
  std::vector<VertexID> candidates;
  VertexID filter_offset = -1;
  VertexID filter_end = 0;

  FilterContext(VertexID offset, VertexID end) : filter_offset(offset), filter_end(end) {}
};

class Filter {
 protected:
  const QueryGraph* query_graph_;
  const Graph* data_graph_;
  const QueryVertexID query_vertex_;
  const std::vector<QueryVertexID> pivot_vertices_;

 private:
  inline bool intersection_not_null_bs(const std::pair<const VertexID*, uint32_t>& set1,
                                       const std::pair<const VertexID*, uint32_t>& set2) {
    if (set1.second <= set2.second) {
      auto lower_bound = set2.first;
      for (uint32_t i = 0; i < set1.second; ++i) {
        auto vid = set1.first[i];
        lower_bound = std::lower_bound(lower_bound, set2.first + set2.second, vid);
        uint32_t index = lower_bound - set2.first;

        if (index >= set2.second) {
          return false;
        }
        if (*lower_bound == vid) {
          return true;
        }
      }
    } else {
      intersection_not_null_bs(set2, set1);
    }
  }

  inline bool intersection_not_null(const std::pair<const VertexID*, uint32_t>& set1,
                                    const std::pair<const VertexID*, uint32_t>& set2) {
    if (std::max(set1.second, set2.second) >= 32) {
      return intersection_not_null_bs(set1, set2);
    }
    if (set1.second <= set2.second) {
      uint32_t set2_index = 0;
      for (uint32_t i = 0; i < set1.second; ++i) {
        auto vid = set1.first[i];
        while (set2_index < set2.second && set2.first[set2_index] < vid) {
          ++set2_index;
        }
        if (set2_index == set2.second) {  // all elements in the rest of set2 are smaller than the rest of set1
          return false;
        }

        if (set2.first[set2_index] == vid) {
          return true;
        }
      }
    } else {
      intersection_not_null(set2, set1);
    }
  }

 public:
  Filter(const QueryGraph* query_graph, const Graph* data_graph, const QueryVertexID query_vertex)
      : query_graph_(query_graph), data_graph_(data_graph), query_vertex_(query_vertex) {}

  Filter(const QueryGraph* query_graph, const Graph* data_graph, const QueryVertexID query_vertex,
         const std::vector<QueryVertexID>& pivot_vertices)
      : query_graph_(query_graph),
        data_graph_(data_graph),
        query_vertex_(query_vertex),
        pivot_vertices_(pivot_vertices) {}

  Filter(const QueryGraph* query_graph, const Graph* data_graph, const QueryVertexID query_vertex,
         const QueryVertexID& pivot_vertex)
      : query_graph_(query_graph),
        data_graph_(data_graph),
        query_vertex_(query_vertex),
        pivot_vertices_({pivot_vertex}) {}

  void filter(const std::vector<std::vector<VertexID>>& candidates, FilterContext* ctx) {
    for (uint32_t i = ctx->offset; i < ctx->end; ++i) {
      VertexID candidate = candidates[query_vertex_][i];
      const auto candidate_nbrs = data_graph_->getOutNeighbors(candidate);
      bool still_candidate = true;
      for (QueryVertexID pivot_vertex : pivot_vertices_) {
        still_candidate &= intersection_not_null(
            candidate_nbrs, std::make_pair(candidates[pivot_vertex].data(), (uint32_t)candidates[pivot_vertex].size()));
      }
      if (still_candidate) {
        ctx->candidates.emplace_back(candidate);
      }
    }
  }
};

}  // namespace circinus

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

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "ops/operator.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

namespace circinus {

struct FilterContext {
  uint32_t offset = -1;
  uint32_t end = 0;

  FilterContext(VertexID offset, VertexID end) : offset(offset), end(end) {}
};

class NeighborhoodFilter : public Operator {
 protected:
  const QueryGraph* query_graph_;
  const QueryVertexID query_vertex_;
  const std::vector<QueryVertexID> pivot_vertices_;
  std::string name_;
  VertexID filter_size_ = 0;

 private:
  inline bool intersection_not_null_bs(const std::pair<const VertexID*, uint32_t>& set1,
                                       const std::pair<const VertexID*, uint32_t>& set2) const {
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
      return intersection_not_null_bs(set2, set1);
    }
  }

  inline bool intersection_not_null(const std::pair<const VertexID*, uint32_t>& set1,
                                    const std::pair<const VertexID*, uint32_t>& set2) const {
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
      return intersection_not_null(set2, set1);
    }
  }

 public:
  NeighborhoodFilter(ExecutionConfig& conf, const QueryGraph* query_graph, const QueryVertexID query_vertex,
                     std::string&& name = "NeighborhoodFilter")
      : Operator(conf.getParallelism()),
        filter_size_(conf.getInputSize()),
        query_graph_(query_graph),
        query_vertex_(query_vertex),
        name_(std::move(name)) {}

  NeighborhoodFilter(ExecutionConfig& conf, const QueryGraph* query_graph, const QueryVertexID query_vertex,
                     const std::vector<QueryVertexID>& pivot_vertices, std::string&& name = "NeighborhoodFilter")
      : Operator(conf.getParallelism()),
        filter_size_(conf.getInputSize()),
        query_graph_(query_graph),
        query_vertex_(query_vertex),
        pivot_vertices_(pivot_vertices),
        name_(std::move(name)) {}

  NeighborhoodFilter(ExecutionConfig& conf, const QueryGraph* query_graph, const QueryVertexID query_vertex,
                     const QueryVertexID& pivot_vertex, std::string&& name = "NeighborhoodFilter")
      : Operator(conf.getParallelism()),
        filter_size_(conf.getInputSize()),
        query_graph_(query_graph),
        query_vertex_(query_vertex),
        pivot_vertices_({pivot_vertex}),
        name_(std::move(name)) {}

  inline uint32_t getParallelism() const { return parallelism_; }
  inline QueryVertexID getQueryVertex() const { return query_vertex_; }
  inline FilterContext initFilterContext(uint32_t task_idx) const {
    DCHECK_LT(task_idx, parallelism_);
    auto chunk_size = filter_size_ / parallelism_;
    CHECK_NE(chunk_size, 0) << "scan_size=" << filter_size_ << ", parallelism=" << parallelism_;
    if (task_idx < filter_size_ % parallelism_) {
      return FilterContext((chunk_size + 1) * task_idx, chunk_size + 1);
    }
    return FilterContext(chunk_size * task_idx + (filter_size_ % parallelism_), chunk_size);
  }

  virtual void filter(const Graph* g, std::vector<std::vector<VertexID>>* candidates, FilterContext* ctx) const {
    for (uint32_t i = ctx->offset; i < ctx->end; ++i) {
      VertexID& candidate = (*candidates)[query_vertex_][i];
      const auto candidate_nbrs = g->getOutNeighbors(candidate);
      bool still_candidate = true;
      for (QueryVertexID pivot_vertex : pivot_vertices_) {
        still_candidate &= intersection_not_null(
            candidate_nbrs,
            std::make_pair((*candidates)[pivot_vertex].data(), (uint32_t)(*candidates)[pivot_vertex].size()));
        if (!still_candidate) {
          candidate = INVALID_VERTEX_ID;
        }
      }
    }
  }
};

}  // namespace circinus

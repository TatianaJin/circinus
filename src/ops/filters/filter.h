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
#include <string>
#include <utility>
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
  // TODO(tatiana): use template forward iterator
  static inline bool intersectionNotNullBS(const std::pair<const VertexID*, uint32_t>& set1,
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
      return intersectionNotNullBS(set2, set1);
    }
    return false;
  }

 public:
  NeighborhoodFilter(ExecutionConfig& conf, const QueryGraph* query_graph, const QueryVertexID query_vertex,
                     std::string&& name = "NeighborhoodFilter")
      : Operator(conf.getParallelism()),
        query_graph_(query_graph),
        query_vertex_(query_vertex),
        name_(std::move(name)),
        filter_size_(conf.getInputSize()) {}

  NeighborhoodFilter(ExecutionConfig& conf, const QueryGraph* query_graph, const QueryVertexID query_vertex,
                     const std::vector<QueryVertexID>& pivot_vertices, std::string&& name = "NeighborhoodFilter")
      : Operator(conf.getParallelism()),
        query_graph_(query_graph),
        query_vertex_(query_vertex),
        pivot_vertices_(pivot_vertices),
        name_(std::move(name)),
        filter_size_(conf.getInputSize()) {}

  NeighborhoodFilter(ExecutionConfig& conf, const QueryGraph* query_graph, const QueryVertexID query_vertex,
                     const QueryVertexID& pivot_vertex, std::string&& name = "NeighborhoodFilter")
      : NeighborhoodFilter(conf, query_graph, query_vertex, std::vector<QueryVertexID>{pivot_vertex}, std::move(name)) {
  }

  static inline bool intersectionNotNull(const std::pair<const VertexID*, uint32_t>& set1,
                                         const std::pair<const VertexID*, uint32_t>& set2) {
    if (std::max(set1.second, set2.second) >= 32) {
      return intersectionNotNullBS(set1, set2);
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
      return intersectionNotNull(set2, set1);
    }
    return false;
  }

  inline void setInputSize(uint64_t filter_size) { filter_size_ = filter_size; }

  inline QueryVertexID getQueryVertex() const { return query_vertex_; }

  inline FilterContext initFilterContext(uint32_t task_idx) const {
    DCHECK_LT(task_idx, parallelism_);
    auto chunk_size = filter_size_ / parallelism_;
    CHECK_NE(chunk_size, 0) << "scan_size=" << filter_size_ << ", parallelism=" << parallelism_;
    if (task_idx < filter_size_ % parallelism_) {
      auto offset = (chunk_size + 1) * task_idx;
      return FilterContext(offset, offset + chunk_size + 1);
    }
    LOG(INFO) << filter_size_;
    auto offset = chunk_size * task_idx + (filter_size_ % parallelism_);
    return FilterContext(offset, offset + chunk_size);
  }

  virtual void filter(const GraphBase* g, std::vector<std::vector<VertexID>>* candidates, FilterContext* ctx) const {
    CHECK_LE(ctx->end, (*candidates)[query_vertex_].size()) << ctx->end << "  " << (*candidates)[query_vertex_].size();
    for (uint32_t i = ctx->offset; i < ctx->end; ++i) {
      VertexID& candidate = (*candidates)[query_vertex_][i];
      const auto candidate_nbrs = g->getOutNeighbors(candidate);
      for (QueryVertexID pivot_vertex : pivot_vertices_) {
        if (!intersectionNotNull(candidate_nbrs, std::make_pair((*candidates)[pivot_vertex].data(),
                                                                (uint32_t)(*candidates)[pivot_vertex].size()))) {
          candidate = INVALID_VERTEX_ID;
          break;
        }
      }
    }
  }
};

}  // namespace circinus

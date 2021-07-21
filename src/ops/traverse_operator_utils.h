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

#include <utility>

#include "graph/bipartite_graph.h"
#include "graph/graph.h"
#include "graph/graph_view.h"
#include "graph/partitioned_graph.h"
#include "graph/types.h"
#include "ops/traverse_operator.h"

namespace circinus {

template <template <typename> typename OpType, typename... Args>
inline TraverseOperator* newTraverseOp(GraphType g_type, Args... args) {
  switch (g_type) {
  case GraphType::Normal:
    return new OpType<Graph>(std::forward<Args>(args)...);
  case GraphType::Partitioned:
    return new OpType<ReorderedPartitionedGraph>(std::forward<Args>(args)...);
  case GraphType::GraphView:
    return new OpType<GraphView<GraphPartitionBase>>(std::forward<Args>(args)...);
  case GraphType::BipartiteGraphView:
    return new OpType<GraphView<BipartiteGraph>>(std::forward<Args>(args)...);
  }
  LOG(FATAL) << "Unknown graph type " << ((uint32_t)g_type);
  return nullptr;
}

template <template <typename, bool> typename OpType, typename... Args>
inline TraverseOperator* newTraverseOp(GraphType g_type, Args... args) {
  switch (g_type) {
  case GraphType::Normal:
    return new OpType<Graph, true>(std::forward<Args>(args)...);
  case GraphType::Partitioned:
    return new OpType<ReorderedPartitionedGraph, true>(std::forward<Args>(args)...);
  case GraphType::GraphView:
    return new OpType<GraphView<GraphPartitionBase>, true>(std::forward<Args>(args)...);
  case GraphType::BipartiteGraphView:
    return new OpType<GraphView<BipartiteGraph>, false>(std::forward<Args>(args)...);
  }
  LOG(FATAL) << "Unknown graph type " << ((uint32_t)g_type);
  return nullptr;
}

template <typename G>
inline constexpr bool sensitive_to_hint = !std::is_same<G, Graph>::value;

inline void removeExceptions(std::vector<VertexID>* set, const unordered_set<VertexID>& except = {}) {
  if (except.empty()) return;
  set->erase(std::remove_if(set->begin(), set->end(), [&except](VertexID v) { return except.count(v); }), set->end());
}

class IntersectionCache {
 private:
  VertexID cache_key_ = std::numeric_limits<VertexID>::max();
  std::vector<VertexID> cache_value_;

 public:
  bool getIntersectionCache(VertexID parent, std::vector<VertexID>* res) {
    if (parent == cache_key_) {
      res->insert(res->end(), cache_value_.begin(), cache_value_.end());
      return true;
    }
    return false;
  }

  void cacheIntersection(VertexID parent, const std::vector<VertexID>& res) {
    cache_key_ = parent;
    cache_value_ = res;
  }
};

class MultiparentIntersectionCache {
 private:
  std::vector<VertexID> cache_keys_;
  std::vector<std::vector<VertexID>> cache_values_;

 public:
  inline void initCacheSize(uint32_t size) {
    cache_keys_.resize(size, std::numeric_limits<VertexID>::max());
    cache_values_.resize(size);
  }

  inline bool intersectionIsNotCached(VertexID key, uint32_t cache_idx) { return cache_keys_[cache_idx] != key; }

  inline const auto& getIntersectionCache(uint32_t cache_idx) const { return cache_values_[cache_idx]; }

  inline auto& resetIntersectionCache(uint32_t cache_idx, VertexID key) {
    CHECK_LT(cache_idx, cache_values_.size());
    cache_values_[cache_idx].clear();
    cache_keys_[cache_idx] = key;
    return cache_values_[cache_idx];
  }

  inline void invalidateCache(uint32_t cache_idx) {
    for (uint32_t i = cache_idx; i < cache_keys_.size(); ++i) {
      cache_keys_[i] = std::numeric_limits<VertexID>::max();
    }
  }
};

}  // namespace circinus

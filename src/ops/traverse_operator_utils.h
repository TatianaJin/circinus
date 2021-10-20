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

#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/intersect.h"
#include "graph/bipartite_graph.h"
#include "graph/graph.h"
#include "graph/graph_view.h"
#include "graph/partitioned_graph.h"
#include "graph/types.h"
#include "graph/vertex_set_view.h"
#include "ops/filters/target_filter.h"
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
[[deprecated]] inline TraverseOperator* newTraverseOp(GraphType g_type, Args... args) {
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

template <template <typename, bool> typename OpType, typename... Args>
inline TraverseOperator* newTraverseOp(GraphType g_type, bool intersect_candidates, Args... args) {
  if (intersect_candidates) {
    switch (g_type) {
    case GraphType::Normal:
      return new OpType<Graph, true>(std::forward<Args>(args)...);
    case GraphType::Partitioned:
      return new OpType<ReorderedPartitionedGraph, true>(std::forward<Args>(args)...);
    case GraphType::GraphView:
      return new OpType<GraphView<GraphPartitionBase>, true>(std::forward<Args>(args)...);
    case GraphType::BipartiteGraphView:
      return new OpType<GraphView<BipartiteGraph>, false>(std::forward<Args>(args)...);  // no need to intersect
    }
  } else {
    switch (g_type) {
    case GraphType::Normal:
      return new OpType<Graph, true>(std::forward<Args>(args)...);  // must intersect to ensure correctness now
    case GraphType::Partitioned:
      return new OpType<ReorderedPartitionedGraph, false>(std::forward<Args>(args)...);
    case GraphType::GraphView:
      return new OpType<GraphView<GraphPartitionBase>, false>(std::forward<Args>(args)...);
    case GraphType::BipartiteGraphView:
      return new OpType<GraphView<BipartiteGraph>, false>(std::forward<Args>(args)...);
    }
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

inline void removeExceptions(const SingleRangeVertexSetView& set, std::vector<VertexID>* res,
                             const unordered_set<VertexID>& except) {
  for (auto vid : set) {
    if (except.count(vid) == 0) res->emplace_back(vid);
  }
}

inline void removeExceptions(const VertexSetView& set, std::vector<VertexID>* res,
                             const unordered_set<VertexID>& except = {}) {
  auto& ranges = set.getRanges();
  for (auto& range : ranges) {
    for (uint32_t i = 0; i < range.second; ++i) {
      auto vid = range.first[i];
      if (except.count(vid) == 0) res->emplace_back(vid);
    }
  }
}

/* dummy functions >>>>> */
template <typename G>
inline void degreeFilter(const VertexSetView& set, uint32_t degree, const GraphView<G>* graph,
                         std::vector<VertexID>* res, const unordered_set<VertexID>& except) {
  removeExceptions(set, res, except);
}

template <typename G>
inline void degreeFilter(std::vector<VertexID>* res, uint32_t degree, const GraphView<G>* graph) {}
/* <<<<< dummy functions */

/** Copy vertices with degree no less than the given degree to res. */
inline void degreeFilter(const SingleRangeVertexSetView& set, uint32_t degree, const GraphBase* graph,
                         std::vector<VertexID>* res, const unordered_set<VertexID>& except) {
  for (auto vid : set) {
    if (graph->getVertexOutDegree(vid) >= degree && except.count(vid) == 0) {
      res->push_back(vid);
    }
  }
}

/** Copy vertices with degree no less than the given degree to res. */
inline void degreeFilter(const VertexSetView& set, uint32_t degree, const GraphBase* graph, std::vector<VertexID>* res,
                         const unordered_set<VertexID>& except) {
  auto& ranges = set.getRanges();
  for (auto& range : ranges) {
    for (uint32_t i = 0; i < range.second; ++i) {
      auto vid = range.first[i];
      if (graph->getVertexOutDegree(vid) >= degree && except.count(vid) == 0) {
        res->push_back(vid);
      }
    }
  }
}

/** Copy vertices with degree no less than the given degree to res. */
inline void degreeFilter(std::vector<VertexID>* res, uint32_t degree, const GraphBase* graph) {
  res->erase(std::remove_if(res->begin(), res->end(),
                            [degree, graph](VertexID v) { return graph->getVertexOutDegree(v) < degree; }),
             res->end());
}

template <typename NeighborSet>
inline void intersectCandidateSetWithProfile(const CandidateSetView& candidates, const NeighborSet& neighbors,
                                             std::vector<VertexID>* targets, const unordered_set<VertexID>& exceptions,
                                             TraverseContext* ctx, const CompressedSubgraphs& input,
                                             TargetFilter* tfilter) {
  removeExceptions(neighbors, targets, exceptions);
  {  // enforce partial order
    if (targets->empty()) return;
    tfilter->filter(targets, input);
    if (targets->empty()) return;
  }
  auto target_size = targets->size();
  intersectInplace(*targets, candidates, targets);
  ctx->candidate_si_diff += target_size - targets->size();
  ctx->updateIntersectInfo(candidates.size() + target_size, targets->size());
}

class TargetBuffer {
 private:
  std::vector<VertexID> current_targets_;  // calculated from current_inputs_[input_index_]
  uint32_t current_target_index_ = 0;
  SingleRangeVertexSetView targets_view_;

 public:
  // for key to key
  inline bool hasTarget() const {
    DCHECK(current_targets_.empty() || current_targets_.size() == targets_view_.size())
        << "forgot to call resetTargetView()?";
    return current_target_index_ < targets_view_.size();
  }

  inline VertexID currentTarget() const { return targets_view_[current_target_index_]; }
  inline void nextTarget() { ++current_target_index_; }

  inline void setTargetView(const SingleRangeVertexSetView& view) {
    targets_view_ = view;
    resetTargets();
  }

  inline void resetTargetView() {
    if (current_targets_.empty()) return;
    targets_view_ = SingleRangeVertexSetView(current_targets_.data(), current_targets_.size());
  }

  /** Should be used in combination with resetTargetView() */
  inline std::vector<VertexID>& resetTargets() {
    current_target_index_ = 0;
    current_targets_.clear();
    targets_view_.clear();
    return current_targets_;
  }
};

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

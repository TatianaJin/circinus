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

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph/types.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/traverse_operator.h"
#include "ops/traverse_operator_utils.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

#ifdef INTERSECTION_CACHE
class ExpandEdgeTraverseContext : public TraverseContext, public IntersectionCache {
#else
class ExpandEdgeTraverseContext : public TraverseContext {
#endif
  unordered_set<VertexID> parent_set_;  // for profile

 public:
  ExpandEdgeTraverseContext(const std::vector<CompressedSubgraphs>* inputs, const void* data_graph,
                            uint32_t input_index, uint32_t input_end_index)
      : TraverseContext(inputs, data_graph, input_index, input_end_index) {}

  template <QueryType profile>
  inline void updateIntersection(uint32_t input_size, uint32_t output_size, VertexID parent) {
    if
      constexpr(isProfileMode(profile)) {
        if
          constexpr(isProfileWithMiniIntersectionMode(profile)) {
            distinct_intersection_count += parent_set_.insert(parent).second;
          }
        total_intersection_input_size += input_size;
        total_intersection_output_size += output_size;
      }
  }

  inline bool hasIntersectionParent(VertexID parent_match) { return parent_set_.insert(parent_match).second; }
};

class ExpandEdgeOperator : public TraverseOperator {
 protected:
  uint32_t parent_index_;  // index of parent query vertex in the compressed subgraphs
  uint32_t target_index_;  // index of target query vertex in the compressed subgraphs
  QueryVertexID parent_id_;
  LabelID parent_label_ = ALL_LABEL;

 public:
  /**
   * Create a new TraverseOperator, which needs to be deleted manually, according to whether the parent and target
   * vertex are compression key or not.
   *
   * @param indices The map from the id of a query vertex to its index in the CompressedSubgraphs
   */
  static TraverseOperator* newExpandEdgeKeyToKeyOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold, SubgraphFilter* filter,
                                                         GraphType graph_type);

  static TraverseOperator* newExpandEdgeKeyToSetOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold, SubgraphFilter* filter,
                                                         GraphType graph_type);

  static TraverseOperator* newExpandEdgeSetToKeyOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold, SubgraphFilter* filter,
                                                         GraphType graph_type);

  ExpandEdgeOperator(uint32_t parent_index, uint32_t target_index, QueryVertexID parent, QueryVertexID target,
                     const std::vector<uint32_t>& same_label_key_indices,
                     const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                     SubgraphFilter* filter)
      : TraverseOperator(target, same_label_key_indices, same_label_set_indices, set_pruning_threshold, filter),
        parent_index_(parent_index),
        target_index_(target_index),
        parent_id_(parent) {}

  virtual ~ExpandEdgeOperator() {}

  inline void setParentLabel(LabelID l) { parent_label_ = l; }

  std::unique_ptr<TraverseContext> initTraverseContext(const std::vector<CompressedSubgraphs>* inputs,
                                                       const void* graph, uint32_t start, uint32_t end,
                                                       QueryType profile) const override {
    auto ret = std::make_unique<ExpandEdgeTraverseContext>(inputs, graph, start, end);
    ret->query_type = profile;
    return ret;
  }

  std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<CandidateSetView>& candidate_sets) override {
    std::vector<std::unique_ptr<BipartiteGraph>> ret;
    ret.reserve(1);
    ret.emplace_back(std::make_unique<BipartiteGraph>(parent_id_, target_vertex_));
    ret.back()->populateGraph(g, candidate_sets);
    return ret;
  }

  std::vector<std::unique_ptr<GraphPartitionBase>> computeGraphPartitions(
      const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const override {
    std::vector<std::unique_ptr<GraphPartitionBase>> ret;
    ret.emplace_back(
        GraphPartitionBase::createGraphPartition(candidate_scopes[parent_id_], candidate_scopes[target_vertex_], g));
    return ret;
  }

  template <typename G, QueryType profile>
  inline void expandFromParent(TraverseContext* ctx, VertexID parent_match,
                               const unordered_set<VertexID>& candidate_set, const unordered_set<VertexID>& exceptions,
                               std::vector<VertexID>* targets) const {
    auto g = (const G*)(ctx->current_data_graph);
#ifdef INTERSECTION_CACHE
    auto dctx = dynamic_cast<ExpandEdgeTraverseContext*>(ctx);
    if (dctx->getIntersectionCache(parent_match, targets)) {
      if
        constexpr(isProfileMode(profile)) { ++ctx->cache_hit; }
    } else {
      auto neighbors = g->getOutNeighborsWithHint(parent_match, target_label_, 0);
      intersect(candidate_set, neighbors, targets);
      dctx->updateIntersection<profile>(candidate_set.size() + neighbors.size(), targets->size(), parent_match);
      dctx->cacheIntersection(parent_match, *targets);
    }
    removeExceptions(targets, exceptions);
#else
    auto neighbors = g->getOutNeighborsWithHint(parent_match, target_label_, 0);
    intersect(candidate_set, neighbors, targets, exceptions);
    ((ExpandEdgeTraverseContext*)ctx)
        ->updateIntersection<profile>(candidate_set.size() + neighbors.size(), targets->size(), parent_match);
#endif
  }

 protected:
  inline void toStringInner(std::stringstream& ss) const {
    ss << ' ' << parent_id_ << ':' << parent_label_ << " -> " << target_vertex_ << ':' << target_label_;
    if (candidates_ != nullptr) ss << " (" << candidates_->size() << ")";
  }
};

}  // namespace circinus

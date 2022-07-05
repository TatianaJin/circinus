#pragma once

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/intersect.h"
#include "algorithms/leapfrog_join.h"
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
#ifdef USE_LFJ
  const CandidateSetView* candidates_;
#else
  unordered_set<VertexID> candidate_set_;
  const unordered_set<VertexID>* candidate_hashmap_;
#endif
  unordered_set<VertexID> parent_set_;  // for profile

 public:
  ExpandEdgeTraverseContext(const CandidateSetView* candidates, const void* graph,
                            std::vector<CompressedSubgraphs>* outputs, QueryType profile, bool hash_candidates,
                            const unordered_set<VertexID>* candidate_hashmap)
      : TraverseContext(graph, outputs, profile) {
#ifdef USE_LFJ
    candidates_ = candidates;
#else
    if (hash_candidates) {
      if (candidate_hashmap != nullptr) {
        candidate_hashmap_ = candidate_hashmap;
      } else {
        candidate_set_ = unordered_set<VertexID>(candidates->begin(), candidates->end());
        candidate_hashmap_ = &candidate_set_;
      }
    }
#endif
  }

  virtual ~ExpandEdgeTraverseContext() {}

  std::unique_ptr<TraverseContext> clone() const override { return std::make_unique<ExpandEdgeTraverseContext>(*this); }

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

#ifdef USE_LFJ
  inline SingleRangeVertexSetView getCandidateSet() const {
    return SingleRangeVertexSetView(candidates_->begin(), candidates_->end());
  }

  inline auto getCandidateSize() const { return candidates_->size(); }
#else
  inline const unordered_set<VertexID>& getCandidateSet() const { return *candidate_hashmap_; }
  inline auto getCandidateSize() const { return candidate_hashmap_->size(); }
#endif
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
                                                         uint64_t set_pruning_threshold,
                                                         std::unique_ptr<SubgraphFilter>&& filter, GraphType graph_type,
                                                         bool intersect_candidates);

  static TraverseOperator* newExpandEdgeKeyToSetOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold,
                                                         std::unique_ptr<SubgraphFilter>&& filter, GraphType graph_type,
                                                         bool intersect_candidates);

  static TraverseOperator* newExpandEdgeSetToKeyOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold,
                                                         std::unique_ptr<SubgraphFilter>&& filter, GraphType graph_type,
                                                         bool intersect_candidates);

  ExpandEdgeOperator(uint32_t parent_index, uint32_t target_index, QueryVertexID parent, QueryVertexID target,
                     const std::vector<uint32_t>& same_label_key_indices,
                     const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                     std::unique_ptr<SubgraphFilter>&& filter)
      : TraverseOperator(target, same_label_key_indices, same_label_set_indices, set_pruning_threshold,
                         std::move(filter)),
        parent_index_(parent_index),
        target_index_(target_index),
        parent_id_(parent) {}

  virtual ~ExpandEdgeOperator() {}

  inline void setParentLabel(LabelID l) { parent_label_ = l; }

  bool extend_vertex() const override { return true; }

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

  template <typename G, QueryType profile, bool intersect_candidates>
  inline void expandFromParent(TraverseContext* ctx, VertexID parent_match, const unordered_set<VertexID>& exceptions,
                               std::vector<VertexID>* targets, const CompressedSubgraphs& input) const {
    auto g = ctx->getDataGraph<G>();
    if
      constexpr(!intersect_candidates) {
        auto neighbors = g->getOutNeighborsWithHint(parent_match, target_label_, 0);
        {  // enforce partial order
          filterTargets(neighbors, input);
          if (neighbors.size() != 0) {
            if (target_degree_ == 1) {
              removeExceptions(neighbors, targets, exceptions);
            } else {
              degreeFilter(neighbors, target_degree_, g, targets, exceptions);
            }
          }
        }
        return;
      }

    auto dctx = dynamic_cast<ExpandEdgeTraverseContext*>(ctx);
#ifdef INTERSECTION_CACHE
    if (dctx->getIntersectionCache(parent_match, targets)) {
      if
        constexpr(isProfileMode(profile)) { ++ctx->cache_hit; }
    } else {
      auto neighbors = g->getOutNeighborsWithHint(parent_match, target_label_, 0);
      intersect(dctx->getCandidateSet(), neighbors, targets);
      dctx->updateIntersection<profile>(dctx->getCandidateSet().size() + neighbors.size(), targets->size(),
                                        parent_match);
      dctx->cacheIntersection(parent_match, *targets);
    }
    removeExceptions(targets, exceptions);
    {  // enforce partial order
      filterTargets(targets, input);
    }
#else
    auto neighbors = g->getOutNeighborsWithHint(parent_match, target_label_, 0);
    if
      constexpr(isProfileCandidateSIEffect(profile)) {
        removeExceptions(neighbors, targets, exceptions);
        {  // enforce partial order
          filterTargets(targets, input);
          if (targets->empty()) return;
        }
        auto target_size = targets->size();
#ifdef USE_LFJ
        intersectInplace(*targets, dctx->getCandidateSet(), targets);
#else
        intersectInplace(targets, dctx->getCandidateSet());
#endif
        ctx->candidate_si_diff += target_size - targets->size();
        ((ExpandEdgeTraverseContext*)ctx)
            ->updateIntersection<profile>(dctx->getCandidateSet().size() + target_size, targets->size(), parent_match);
        return;
      }
#ifdef USE_LFJ
    leapfrogJoin(dctx->getCandidateSet(), neighbors, targets, exceptions);
#else
    intersect(dctx->getCandidateSet(), neighbors, targets, exceptions);
#endif
    if
      constexpr(isProfileMode(profile)) {
        ((ExpandEdgeTraverseContext*)ctx)
            ->updateIntersection<profile>(dctx->getCandidateSet().size() + neighbors.size(), targets->size(),
                                          parent_match);
      }
    {  // enforce partial order
      filterTargets(targets, input);
    }
#endif
  }

 protected:
  inline void toStringInner(std::stringstream& ss) const {
    ss << ' ' << parent_id_ << ':' << parent_label_ << " -> " << target_vertex_ << ':' << target_label_;
  }
};

}  // namespace circinus

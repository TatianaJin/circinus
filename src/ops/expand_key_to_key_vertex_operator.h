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
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/expand_vertex_traverse_context.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandKeyToKeyVertexTraverseContext : public ExpandVertexTraverseContext, public TargetBuffer {
 public:
  ExpandKeyToKeyVertexTraverseContext(const CandidateSetView* candidates, const void* graph,
                                      std::vector<CompressedSubgraphs>* outputs, QueryType profile,
                                      uint32_t parent_size)
      : ExpandVertexTraverseContext(candidates, graph, outputs, profile, parent_size) {}

  std::unique_ptr<TraverseContext> clone() const override {
    return std::make_unique<ExpandKeyToKeyVertexTraverseContext>(*this);
  }
};

template <typename G, bool intersect_candidates>
class ExpandKeyToKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandKeyToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               SubgraphFilter* filter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, filter) {
    CHECK_GT(parents.size(), 1);
  }

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, (ExpandKeyToKeyVertexTraverseContext*)ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    if (ctx->getQueryType() == QueryType::Profile) {
      return expandInner<QueryType::Profile>(batch_size, (ExpandKeyToKeyVertexTraverseContext*)ctx);
    }
    if (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect) {
      return expandInner<QueryType::ProfileCandidateSIEffect>(batch_size, (ExpandKeyToKeyVertexTraverseContext*)ctx);
    }
    CHECK(ctx->getQueryType() == QueryType::ProfileWithMiniIntersection) << "Unknown query type "
                                                                         << (uint32_t)ctx->getQueryType();
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, (ExpandKeyToKeyVertexTraverseContext*)ctx);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToKeyVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first + 1, input_key_size.second + 1};
  }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    auto ret =
        std::make_unique<ExpandKeyToKeyVertexTraverseContext>(candidates, graph, outputs, profile, parents_.size());
#ifdef INTERSECTION_CACHE
    ret->initCacheSize(parents_.size());
#endif
    return ret;
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(uint32_t batch_size, ExpandKeyToKeyVertexTraverseContext* ctx) const {
    auto data_graph = ctx->getDataGraph<G>();
    uint32_t output_num = 0;
    while (true) {
      if (ctx->hasTarget()) {
        auto& input = ctx->getPreviousInput();
        while (ctx->hasTarget()) {
#ifdef USE_FILTER
          auto output =
              ctx->newOutput(input, ctx->currentTarget(), same_label_set_indices_, set_pruning_threshold_, false);
          ctx->nextTarget();
          if (output == nullptr) continue;
          if (filter(*output)) {
            ctx->popOutput();
            continue;
          }
#else
          auto output = ctx->newOutput(input, ctx->currentTarget(), same_label_set_indices_, set_pruning_threshold_);
          ctx->nextTarget();
          if (output == nullptr) continue;
#endif
          if (++output_num == batch_size) {
            return output_num;
          }
        }
      }
      if (!ctx->hasNextInput()) {  // return if all inputs in the current batch are consumed
        return output_num;
      }

      // consume the next input
      auto& new_keys = ctx->resetTargets();
      const auto& input = ctx->getCurrentInput();
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);

#ifdef INTERSECTION_CACHE
      uint32_t i = 0;
      for (; i < parents_.size(); ++i) {
        uint32_t key_vid = input.getKeyVal(query_vertex_indices_.at(parents_[i]));
        if (ctx->intersectionIsNotCached(key_vid, i)) {
          break;
        }
      }
      if (i != 0) {
        if
          constexpr(isProfileMode(profile)) ctx->cache_hit += i;
      } else {
        uint32_t key_vid = input.getKeyVal(query_vertex_indices_.at(parents_[0]));
        auto neighbors = data_graph->getOutNeighborsWithHint(key_vid, target_label_, 0);
        if (intersect_candidates) {
          auto& intersection = ctx->resetIntersectionCache(0, key_vid);
          intersect(*ctx->getCandidateSet(), neighbors, &intersection);
          i = 1;
          if
            constexpr(isProfileMode(profile)) {
              ctx->updateIntersectInfo(ctx->getCandidateSet()->size() + neighbors.size(), intersection.size());
            }
        } else {
          uint32_t second_key_vid = input.getKeyVal(query_vertex_indices_.at(parents_[1]));
          auto second_parent_neighbors = data_graph->getOutNeighborsWithHint(second_key_vid, target_label_, 1);
          ctx->resetIntersectionCache(0, key_vid);
          auto& intersection = ctx->resetIntersectionCache(1, second_key_vid);
          intersect(neighbors, second_parent_neighbors, &intersection);
          i = 2;
          if
            constexpr(isProfileMode(profile)) {
              ctx->updateIntersectInfo(neighbors.size() + second_parent_neighbors.size(), intersection.size());
            }
        }
      }
      if (ctx->getIntersectionCache(i - 1).empty() && (intersect_candidates || i != 1)) {
        ctx->invalidateCache(i);
        goto handle_output;
      }
      for (; i < parents_.size(); ++i) {
        uint32_t key_vid = input.getKeyVal(query_vertex_indices_.at(parents_[i]));
        auto neighbors = data_graph->getOutNeighborsWithHint(key_vid, target_label_, i);
        const auto& last = ctx->getIntersectionCache(i - 1);
        auto& current = ctx->resetIntersectionCache(i, key_vid);
        intersect(last, neighbors, &current);
        if
          constexpr(isProfileMode(profile)) {
            ctx->updateIntersectInfo(last.size() + neighbors.size(), current.size());
          }
        if (current.empty()) {
          ctx->invalidateCache(i + 1);
          break;
        }
      }
      if (i == parents_.size()) {
        new_keys = ctx->getIntersectionCache(i - 1);
        removeExceptions(&new_keys, exceptions);
      }
    handle_output:
#else
      expandFromParents<G, profile, intersect_candidates>(input, data_graph, ctx, parent_indices_, exceptions,
                                                          &new_keys);
#endif

      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
          if
            constexpr(isProfileWithMiniIntersectionMode(profile)) {
              // consider reuse of partial intersection results at each parent
              std::vector<VertexID> parent_tuple(parents_.size());
              for (uint32_t j = 0; j < parents_.size(); ++j) {
                uint32_t key_vid = input.getKeyVal(parent_indices_[j]);
                parent_tuple[j] = key_vid;
                if (j != 0 || intersect_candidates) {
                  ctx->updateDistinctSICount(j, parent_tuple, j);
                }
              }
            }
        }
      ctx->nextInput();
    }
    return output_num;
  }
};

}  // namespace circinus

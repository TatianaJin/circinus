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

template <typename G, bool intersect_candidates>
class ExpandKeyToSetVertexOperator : public ExpandVertexOperator {
 private:
  bool have_existing_target_set_ = false;
  std::vector<std::pair<uint32_t, uint32_t>> uncovered_parent_indices_;  // parent key index, parent index

 public:
  ExpandKeyToSetVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               std::unique_ptr<SubgraphFilter>&& sfilter = nullptr,
                               bool have_existing_target_set = false)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, std::move(sfilter)),
        have_existing_target_set_(have_existing_target_set) {}

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    if (ctx->getQueryType() == QueryType::Profile) return expandInner<QueryType::Profile>(batch_size, ctx);
    if (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect)
      return expandInner<QueryType::ProfileCandidateSIEffect>(batch_size, ctx);
    CHECK(ctx->getQueryType() == QueryType::ProfileWithMiniIntersection) << "Unknown query type "
                                                                         << (uint32_t)ctx->getQueryType();
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, ctx);
  }

  inline void setHaveExistingTargetSet() { have_existing_target_set_ = true; }

  bool extend_vertex() const override { return have_existing_target_set_ == false; }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToSetVertexOperator " << (have_existing_target_set_ ? "with existing target set " : "");
    toStringInner(ss);
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first, input_key_size.second + (have_existing_target_set_ == false)};
  }

  void reuseSetForTarget(uint32_t set_index, const std::vector<QueryVertexID>& uncovered_parents) override {
    uncovered_parent_indices_.reserve(uncovered_parents.size());
    unordered_map<QueryVertexID, uint32_t> parent_index;
    for (auto parent : parents_) {
      parent_index.emplace(parent, parent_index.size());
    }
    for (auto parent : uncovered_parents) {
      uncovered_parent_indices_.emplace_back(query_vertex_indices_.at(parent), parent_index.at(parent));
    }
    reusable_set_index_ = set_index;
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(uint32_t batch_size, TraverseContext* base_ctx) const {
    auto ctx = (ExpandVertexTraverseContext*)base_ctx;
    uint32_t output_num = 0;
    // TODO(tatiana): tidy up
    if (have_existing_target_set_) {
      uint32_t target_vertex_index = query_vertex_indices_.at(target_vertex_);
      for (; output_num < batch_size && ctx->hasNextInput(); ctx->nextInput()) {
        const auto& input = ctx->getCurrentInput();
        std::vector<VertexID> target_set;

        using NeighborSet = typename G::NeighborSet;
        std::vector<NeighborSet> neighbor_sets;
        neighbor_sets.reserve(parent_indices_.size());

        auto g = ctx->getDataGraph<G>();
        for (uint32_t i = 0; i < parent_indices_.size(); ++i) {
          DCHECK_LT(parent_indices_[i], input.getNumKeys());
          uint32_t key_vid = input.getKeyVal(parent_indices_[i]);
          neighbor_sets.push_back(g->getOutNeighborsWithHint(key_vid, target_label_, i));
        }
        std::sort(neighbor_sets.begin(), neighbor_sets.end(),
                  [](const NeighborSet& a, const NeighborSet& b) { return a.size() < b.size(); });
        if (neighbor_sets.front().size() == 0) continue;

        intersect(*(input.getSet(target_vertex_index)), neighbor_sets.front(), &target_set);

        uint32_t si_input_size = input.getSet(target_vertex_index)->size() + neighbor_sets[0].size();
        (void)si_input_size;
        if
          constexpr(isProfileMode(profile)) { ctx->updateIntersectInfo(si_input_size, target_set.size()); }

        {  // enforce partial order
          filterTargets(&target_set, input);
          if (target_set.empty()) {
            continue;
          }
        }

        for (uint32_t i = 1; i < parent_indices_.size(); ++i) {
          auto si_input_size = target_set.size() + neighbor_sets[i].size();
          (void)si_input_size;
          intersectInplace(target_set, neighbor_sets[i], &target_set);
          if
            constexpr(isProfileMode(profile)) { ctx->updateIntersectInfo(si_input_size, target_set.size()); }

          if (target_set.empty()) {
            break;
          }
        }

        if (!target_set.empty()) {
          CompressedSubgraphs& output = ctx->copyOutput(input);
#ifdef USE_FILTER
          output.UpdateSets(target_vertex_index, newVertexSet(target_set));
          if (filter(output)) {
            ctx->popOutput();
            continue;
          }
#else
          if (target_set.size() == 1) {
            unordered_set<uint32_t> pruning_set_indices(same_label_set_indices_.begin(), same_label_set_indices_.end());
            if (output.pruneExistingSets(target_set.front(), pruning_set_indices, set_pruning_threshold_)) {
              ctx->popOutput();
              continue;
            }
          }
          output.UpdateSets(target_vertex_index, newVertexSet(target_set));
#endif
          ++output_num;
        }
      }
      return output_num;
    }

    auto data_graph = ctx->getDataGraph<G>();
    for (; output_num < batch_size && ctx->hasNextInput(); ctx->nextInput()) {
      const auto& input = ctx->getCurrentInput();

      VertexSet target_set;
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      if (canReuseSet()) {
        reuseAndExpand<G, profile, intersect_candidates>(input, data_graph, ctx, uncovered_parent_indices_, exceptions,
                                                         target_set);
        if (target_set->empty()) continue;
      } else {
        std::vector<VertexID> new_set;
        expandFromParents<G, profile, intersect_candidates>(input, data_graph, ctx, parent_indices_, exceptions,
                                                            &new_set);
        if
          constexpr(isProfileMode(profile)) {
            ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
            // consider reuse of partial intersection results at each parent
            if (isProfileWithMiniIntersectionMode(profile)) {
              std::vector<VertexID> parent_tuple(parents_.size());
              for (uint32_t j = 0; j < parent_indices_.size(); ++j) {
                uint32_t key_vid = input.getKeyVal(parent_indices_[j]);
                parent_tuple[j] = key_vid;
                if (j != 0 || intersect_candidates) {
                  ctx->updateDistinctSICount(j, parent_tuple, j);
                }
              }
            }
          }
        if (new_set.empty()) continue;
        target_set = newVertexSet(std::move(new_set));
      }
#ifdef USE_FILTER
      auto output = ctx->newOutput(input, std::move(target_set));
      if (filter(*output)) {
        ctx->popOutput();
        continue;
      }
#else
      auto output = ctx->newOutput(input, std::move(target_set), same_label_set_indices_, set_pruning_threshold_);
      if (output == nullptr) {
        continue;
      }
#endif
      ++output_num;
    }
    return output_num;
  }
};

}  // namespace circinus

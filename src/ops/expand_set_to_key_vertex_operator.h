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
#include <queue>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"
#include "ops/expand_vertex_traverse_context.h"
#include "ops/traverse_operator_utils.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandSetToKeyVertexTraverseContext : public ExpandVertexTraverseContext {
 private:
  unordered_set<VertexID> exceptions_;
  bool from_candidate_ = false;  // init to be false since candidate_iter_ is not initialized before first use
  // for fromCandidateStrategy
  CandidateSetView::ConstIterator candidate_iter_;
  // for fromSetNeighborStrategy
  uint32_t min_parent_idx_ = 0;
  uint32_t parent_match_idx_ = 0;
  const std::vector<VertexID>* parent_matches_ = nullptr;  // not owned
  unordered_set<VertexID> visited_;
  std::queue<VertexID> extensions_;

 public:
  ExpandSetToKeyVertexTraverseContext(const CandidateSetView* candidates, const void* graph,
                                      std::vector<CompressedSubgraphs>* outputs, QueryType profile,
                                      uint32_t parent_size)
      : ExpandVertexTraverseContext(candidates, graph, outputs, profile, parent_size) {}

  std::unique_ptr<TraverseContext> clone() const override {
    return std::make_unique<ExpandSetToKeyVertexTraverseContext>(*this);
  }

  bool fromCandidate() const { return from_candidate_; }

  inline const auto& getExceptions() const { return exceptions_; }
  inline auto& resetExceptions() {
    exceptions_.clear();
    return exceptions_;
  }

  /* for fromCandidateStrategy */
  inline void initCandidateIter() {
    candidate_iter_ = candidates_->begin();
    from_candidate_ = true;
  }
  inline auto& getCandidateIter() { return candidate_iter_; }

  /* for fromSetNeighborStrategy */
  inline auto getParentMatches() const { return parent_matches_; }
  inline auto getMinParentIdx() const { return min_parent_idx_; }
  inline auto& getParentMatchIdx() { return parent_match_idx_; }
  inline auto& getExtensions() { return extensions_; }
  inline bool visitOnce(VertexID target) { return visited_.insert(target).second; }
  inline void initParentSet(const std::vector<VertexID>* parent_matches, uint32_t min_parent_idx) {
    parent_matches_ = parent_matches;
    min_parent_idx_ = min_parent_idx;
    from_candidate_ = false;
  }
  inline void clearParentMatches() {
    parent_match_idx_ = 0;
    parent_matches_ = nullptr;
    visited_.clear();
  }
};

template <typename G, bool intersect_candidates>
class ExpandSetToKeyVertexOperator : public ExpandVertexOperator {
  std::vector<LabelID> parent_labels_;

 public:
  ExpandSetToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               std::unique_ptr<SubgraphFilter>&& subgraph_filter = nullptr,
                               std::vector<LabelID>&& parent_labels = {})
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, std::move(subgraph_filter)),
        parent_labels_(std::move(parent_labels)) {
    if (parent_labels_.empty()) parent_labels_.resize(parents.size(), ALL_LABEL);
  }

  inline void setParentLabels(std::vector<LabelID>&& parent_labels) { parent_labels_ = std::move(parent_labels); }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    return std::make_unique<ExpandSetToKeyVertexTraverseContext>(candidates, graph, outputs, profile, parents_.size());
  }

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

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandSetToKeyVertexOperator";
    {
      uint32_t idx = 0;
      for (auto parent : parents_) {
        DCHECK_EQ(query_vertex_indices_.count(parent), 1);
        ss << ' ' << parent << ':' << parent_labels_[idx++];
      }
      DCHECK_EQ(query_vertex_indices_.count(target_vertex_), 1);
      ss << " -> " << target_vertex_ << ':' << target_label_;
    }
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first + 1, input_key_size.second + 1};
  }

 protected:
  template <typename ForwardIter>
  inline bool isInCandidates(VertexID key, ForwardIter& begin, ForwardIter end, const G* graph) const {
    if
      constexpr(intersect_candidates) {
        begin = std::lower_bound(begin, end, key);
        return begin != end && *begin == key;
      }
    if
      constexpr(std::is_base_of_v<GraphBase, G>) {
        return ((const GraphBase*)graph)->getVertexOutDegree(key) >= target_degree_;
      }
    return true;
  }

  std::tuple<uint32_t, uint32_t, uint32_t> getMinimumParent(TraverseContext* ctx) const {
    uint32_t parent = 0, size = 0xFFFFFFFF, parent_idx = 0;
    auto& input = ctx->getCurrentInput();
    uint32_t idx = 0;
    for (auto par : parents_) {
      auto current_set = input.getSet(query_vertex_indices_.at(par));
      uint32_t current_size = 0;
      for (auto set_vertex_id : *current_set) {
        current_size += ctx->getDataGraph<G>()->getVertexOutDegreeWithHint(set_vertex_id, target_label_, idx);
      }
      if (current_size < size) {
        size = current_size;
        parent = par;
        parent_idx = idx;
      }
      ++idx;
    }
    return std::make_tuple(size / parents_.size(), parent, parent_idx);
  }

  template <QueryType profile>
  inline uint32_t expandInner(uint32_t batch_size, TraverseContext* base_ctx) const {
    auto ctx = (ExpandSetToKeyVertexTraverseContext*)base_ctx;
    uint32_t output_num = 0;
    decltype(std::declval<G>().getInNeighborsWithHint(0, 0, 0)) buffer;

    // check remaining targets for the last input
    if (ctx->fromCandidate()) {
      if (ctx->getCandidateIter() != ctx->getCandidateSet()->end())
        output_num = fromCandidateStrategy<profile>(ctx, ctx->getPreviousInput(), batch_size, buffer);
    } else if (ctx->getParentMatches() != nullptr) {
      output_num = fromSetNeighborStrategy<profile>(ctx, ctx->getPreviousInput(), batch_size, buffer);
    }

    if (output_num == batch_size) return output_num;

    // compute targets for the current input and enumerate
    while (ctx->hasNextInput() && output_num < batch_size) {
      auto[min_parent_set_size, min_parent_vertex, min_parent_idx] = getMinimumParent(ctx);
      if
        constexpr(isProfileWithMiniIntersectionMode(profile)) {
          updateDistinctSICount((ExpandVertexTraverseContext*)ctx);
        }
      auto& input = ctx->getCurrentInput();
      input.getExceptions(ctx->resetExceptions(), same_label_key_indices_, same_label_set_indices_);
      if (ctx->getCandidateSet() == nullptr || min_parent_set_size < ctx->getCandidateSet()->size()) {
        ctx->initParentSet(input.getSet(query_vertex_indices_.at(min_parent_vertex)).get(), min_parent_idx);
        output_num += fromSetNeighborStrategy<profile>(ctx, input, batch_size - output_num, buffer);
      } else {
        ctx->initCandidateIter();
        output_num += fromCandidateStrategy<profile>(ctx, input, batch_size - output_num, buffer);
      }
      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
        }
      ctx->nextInput();
    }
    return output_num;
  }

  /** @returns True if target is valid. */
  template <QueryType profile, typename NeighborSet>
  inline bool validateTarget(ExpandSetToKeyVertexTraverseContext* ctx, VertexID target,
                             const CompressedSubgraphs& input, uint32_t source_parent_idx,
                             NeighborSet& key_out_neighbors) const {
    auto new_output = ctx->newOutput(input, target, same_label_set_indices_, set_pruning_threshold_);
    if (new_output == nullptr) {
      return false;
    }
#ifndef USE_FILTER
    // TODO(tatiana): `ExpandSetToKey` requires different groups of same-label indices for the parent sets
    // for active pruning, should use same-label set indices >>>
    unordered_set<uint32_t> set_indices;
    for (uint32_t i = 0; i < new_output->getNumSets(); ++i) {
      set_indices.insert(i);
    }  // <<< for active pruning, should use same-label set indices
#endif
    // TODO(by) hash key_out_neighbors?
    if
      constexpr(!sensitive_to_hint<G>) key_out_neighbors =
          ctx->getDataGraph<G>()->getInNeighborsWithHint(target, ALL_LABEL, source_parent_idx);
    uint32_t parent_idx = 0;
    for (uint32_t set_vid : parents_) {
      std::vector<VertexID> new_set;
      uint32_t id = query_vertex_indices_.at(set_vid);
      if
        constexpr(sensitive_to_hint<G>) {  // if graph provides different neighbors for hints, select neighbors from
                                           // the right graph part
          if
            constexpr(std::is_same_v<G, ReorderedPartitionedGraph>) ctx->getDataGraph<G>()->getInNeighborsWithHint(
                target, parent_labels_[parent_idx], parent_idx, key_out_neighbors);
          else
            key_out_neighbors =
                ctx->getDataGraph<G>()->getInNeighborsWithHint(target, parent_labels_[parent_idx], parent_idx);
        }
      intersect(*input.getSet(id), key_out_neighbors, &new_set);  // No need for exceptions
      if
        constexpr(isProfileMode(profile)) {
          ctx->updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.size(), new_set.size());
        }
      if (new_set.empty()) {
        ctx->popOutput();
        return false;
      }
#ifdef USE_FILTER
      new_output->UpdateSets(id, newVertexSet(new_set));
      if (filter(*new_output)) {
        ctx->popOutput();
        return false;
      }
#else
      if (new_set.size() == 1) {  // actively prune existing sets
        auto pruning_set_indices = set_indices;
        pruning_set_indices.erase(id);
        if (new_output->pruneExistingSets(new_set.front(), pruning_set_indices, ~0u)) {
          ctx->popOutput();
          return false;
        }
      }
      new_output->UpdateSets(id, newVertexSet(new_set));
#endif
      ++parent_idx;
    }

    return true;
  }

  template <QueryType profile, typename NeighborSet>
  uint32_t fromCandidateStrategy(ExpandSetToKeyVertexTraverseContext* ctx, const CompressedSubgraphs& input,
                                 uint32_t batch_size, NeighborSet& buffer) const {
    uint32_t output_num = 0;
    auto& iter = ctx->getCandidateIter();
    while (iter != ctx->getCandidateSet()->end()) {
      auto vid = *iter;
      ++iter;
      if (ctx->getExceptions().count(vid)) {
        continue;
      }
      if (filterTarget(vid, input)) {  // enforce partial order
        continue;
      }
      output_num += validateTarget<profile>(ctx, vid, input, 0, buffer);
      if (output_num == batch_size) {
        break;
      }
    }
    return output_num;
  }

  template <QueryType profile, typename NeighborSet>
  uint32_t fromSetNeighborStrategy(ExpandSetToKeyVertexTraverseContext* ctx, const CompressedSubgraphs& input,
                                   uint32_t batch_size, NeighborSet& buffer) const {
    uint32_t output_num = 0;
    auto min_parent_idx = ctx->getMinParentIdx();
    auto& extensions = ctx->getExtensions();

    // enumerate remaining extensions from last parent_match_idx
    while (!extensions.empty()) {
      auto key_vertex_id = extensions.front();
      extensions.pop();
      output_num += validateTarget<profile>(ctx, key_vertex_id, input, min_parent_idx, buffer);
      if (output_num == batch_size) return output_num;
    }

    auto& parent_match_idx = ctx->getParentMatchIdx();
    auto& parent_matches = *ctx->getParentMatches();
    auto candidate_end = ctx->getCandidateSet()->end();
    auto graph = ctx->getDataGraph<G>();
    while (parent_match_idx < parent_matches.size()) {
      // find extensions and enumerate
      auto vid = parent_matches[parent_match_idx++];
      auto out_neighbors = graph->getOutNeighborsWithHint(vid, target_label_, min_parent_idx);
      auto begin = ctx->getCandidateSet()->begin();
      if (std::is_same_v<decltype(out_neighbors), VertexSetView>) {
        auto& ranges = out_neighbors.getRanges();
        for (uint32_t range_i = 0; range_i < ranges.size(); ++range_i) {
          for (size_t idx = 0; idx < ranges[range_i].second; ++idx) {
            VertexID key_vertex_id = ranges[range_i].first[idx];

            if (ctx->getExceptions().count(key_vertex_id) == 0 && ctx->visitOnce(key_vertex_id) &&
                !filterTarget(key_vertex_id, input)) {
              if (!isInCandidates(key_vertex_id, begin, candidate_end, graph)) {
                if
                  constexpr(isProfileCandidateSIEffect(profile)) ctx->candidate_si_diff += 1;
                continue;
              }
              output_num += validateTarget<profile>(ctx, key_vertex_id, input, min_parent_idx, buffer);

              if (output_num == batch_size) {
                // store the remaining extensions for next batch output
                for (++idx; idx < ranges[range_i].second; ++idx) {
                  auto vid = ranges[range_i].first[idx];

                  if (ctx->getExceptions().count(vid) == 0 && ctx->visitOnce(vid) && !filterTarget(vid, input)) {
                    if (isInCandidates(vid, begin, candidate_end, graph)) {
                      extensions.push(vid);
                    } else {
                      if
                        constexpr(isProfileCandidateSIEffect(profile)) ctx->candidate_si_diff += 1;
                    }
                  }
                }
                for (++range_i; range_i < ranges.size(); ++range_i) {
                  for (size_t idx = 0; idx < ranges[range_i].second; ++idx) {
                    auto vid = ranges[range_i].first[idx];
                    if (ctx->getExceptions().count(vid) == 0 && ctx->visitOnce(vid) && !filterTarget(vid, input)) {
                      if (isInCandidates(vid, begin, candidate_end, graph)) {
                        extensions.push(vid);
                      } else {
                        if
                          constexpr(isProfileCandidateSIEffect(profile)) ctx->candidate_si_diff += 1;
                      }
                    }
                  }
                }
                return output_num;
              }
            }
          }
        }
      } else {
        for (auto iter = out_neighbors.begin(); iter != out_neighbors.end(); ++iter) {
          VertexID key_vertex_id = *iter;
          if (ctx->getExceptions().count(key_vertex_id) == 0 && ctx->visitOnce(key_vertex_id) &&
              !filterTarget(key_vertex_id, input)) {
            if (!isInCandidates(key_vertex_id, begin, candidate_end, graph)) {
              if
                constexpr(isProfileCandidateSIEffect(profile)) ctx->candidate_si_diff += 1;
              continue;
            }
            output_num += validateTarget<profile>(ctx, key_vertex_id, input, min_parent_idx, buffer);
            if (output_num == batch_size) {
              // store the remaining extensions for next batch output
              for (++iter; iter != out_neighbors.end(); ++iter) {
                auto vid = *iter;
                if (ctx->getExceptions().count(vid) == 0 && ctx->visitOnce(vid) && !filterTarget(vid, input)) {
                  if (isInCandidates(vid, begin, candidate_end, graph)) {
                    extensions.push(vid);
                  } else {
                    if
                      constexpr(isProfileCandidateSIEffect(profile)) ctx->candidate_si_diff += 1;
                  }
                }
              }
              return output_num;
            }
          }
        }
      }
    }
    ctx->clearParentMatches();
    return output_num;
  }

  /** Calculate the ideal si count as if we expand parent 1, 2, ... n for n = parents_.size() in normal backtracing
   * implementation. */
  void updateDistinctSICount(ExpandVertexTraverseContext* ctx) const {
    auto& input = ctx->getCurrentInput();
    std::vector<std::vector<VertexID>*> parent_set_ptrs;
    parent_set_ptrs.reserve(parents_.size());
    for (auto parent : parents_) {
      parent_set_ptrs.push_back(input.getSet(query_vertex_indices_.at(parent)).get());
    }
    uint32_t depth = 0, last_depth = parents_.size() - 1;
    std::vector<uint32_t> set_index(parents_.size(), 0);
    std::vector<VertexID> parent_tuple(parents_.size());
    unordered_set<VertexID> prefix_set;
    while (true) {
      while (set_index[depth] < parent_set_ptrs[depth]->size()) {
        auto parent_vid = (*parent_set_ptrs[depth])[set_index[depth]];
        if (prefix_set.count(parent_vid)) {  // the parent tuples with current prefix will be pruned
          ++set_index[depth];
          continue;
        }
        parent_tuple[depth] = parent_vid;
        if (depth != 0 || intersect_candidates) {
          ctx->updateDistinctSICount(depth, parent_tuple, depth);
        }
        if (depth == last_depth) {
          ++set_index[depth];
        } else {
          prefix_set.insert(parent_vid);
          ++depth;
          set_index[depth] = 0;
        }
      }
      if (depth == 0) {
        break;
      }
      --depth;
      prefix_set.erase((*parent_set_ptrs[depth])[set_index[depth]]);
      ++set_index[depth];
    }
  }
};

}  // namespace circinus

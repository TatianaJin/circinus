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
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"
#include "ops/expand_vertex_traverse_context.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

template <typename G>
class ExpandSetToKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandSetToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               SubgraphFilter* subgraph_filter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, subgraph_filter) {}

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, uint32_t query_type, TraverseContext* ctx) const override {
    if (query_type == 1) {
      return expandInner<QueryType::Profile>(batch_size, ctx);
    }
    CHECK_EQ(query_type, 2) << "unknown query type " << query_type;
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, ctx);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandSetToKeyVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

 protected:
  bool isInCandidates(VertexID key) {
    auto lb = std::lower_bound(candidates_->begin(), candidates_->end(), key);
    return lb != candidates_->end() && *lb == key;
  }

  std::tuple<uint32_t, uint32_t, uint32_t> getMinimumParent(TraverseContext* ctx) const {
    uint32_t parent = 0, size = 0xFFFFFFFF, parent_idx = 0;
    auto& input = ctx->getCurrentInput();
    uint32_t idx = 0;
    for (auto par : parents_) {
      auto current_set = input.getSet(query_vertex_indices_[par]);
      uint32_t current_size = 0;
      for (auto set_vertex_id : *current_set) {
        current_size += ((G*)(ctx->current_data_graph))->getVertexOutDegreeWithHint(set_vertex_id, ALL_LABEL, idx);
      }
      if (current_size < size) {
        size = current_size;
        parent = par;
        parent_idx = idx;
      }
      ++idx;
    }
    return std::make_tuple(size, parent, parent_idx);
  }

  // TODO(tatiana): see if hard limit on output size is needed
  template <QueryType profile>
  inline uint32_t expandInner(uint32_t batch_size, TraverseContext* ctx) const {
    uint32_t output_num = 0;
    while (ctx->hasNextInput()) {
      auto[min_parent_set_size, min_parent_vertex, min_parent_idx] = getMinimumParent();
      if
        constexpr(isProfileWithMiniIntersectionMode(profile)) { updateDistinctSICount(); }
      if (min_parent_set_size < candidates_->size()) {
        DLOG(INFO) << "fromSetNeighborStrategy";
        output_num += fromSetNeighborStrategy<profile>(min_parent_vertex, min_parent_idx, ctx);
      } else {
        DLOG(INFO) << "fromCandidateStrategy";
        output_num += fromCandidateStrategy<profile>(ctx);
      }
      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
        }
      ctx->nextInput();
      if (output_num >= batch_size) {
        break;
      }
    }
    return output_num;
  }

  template <QueryType profile>
  uint32_t fromCandidateStrategy(TraverseContext* ctx) const {
    auto& input = ctx->getCurrentInput();
    uint32_t output_num = 0;
    auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
    for (VertexID key_vertex_id : *candidates_) {
      if (exceptions.count(key_vertex_id)) {
        continue;
      }
      auto key_out_neighbors = ((G*)(ctx->current_data_graph))->getOutNeighborsWithHint(key_vertex_id, ALL_LABEL, 0);
      // TODO(by) hash key_out_neighbors
      CompressedSubgraphs new_output(input, key_vertex_id, same_label_set_indices_, set_pruning_threshold_);
      if (new_output.empty()) {
        continue;
      }
      bool add = true;
#ifndef USE_FILTER
      // TODO(tatiana): `ExpandSetToKey` requires different groups of same-label indices for the parent sets
      // for active pruning, should use same-label set indices >>>
      unordered_set<uint32_t> set_indices;
      for (uint32_t i = 0; i < new_output.getNumSets(); ++i) {
        set_indices.insert(i);
      }  // <<< for active pruning, should use same-label set indices
#endif
      uint32_t parent_idx = 0;
      for (uint32_t set_vid : parents_) {
        std::vector<VertexID> new_set;
        uint32_t id = query_vertex_indices_[set_vid];
        if
          constexpr(
              !std::is_same<G, Graph>::value) {  // if using graph view, select neighbors from the right graph part
            key_out_neighbors =
                ((G*)(ctx->current_data_graph))->getOutNeighborsWithHint(key_vertex_id, ALL_LABEL, parent_idx);
          }
        intersect(*input.getSet(id), key_out_neighbors, &new_set);  // No need for exceptions
        if
          constexpr(isProfileMode(profile)) {
            updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.size(), new_set.size());
          }
        if (new_set.empty()) {
          add = false;
          break;
        }
#ifdef USE_FILTER
        new_output.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
        if (filter(new_output)) {
          add = false;
          break;
        }
#else
        if (new_set.size() == 1) {  // actively prune existing sets
          set_indices.erase(id);
          if (new_output.pruneExistingSets(new_set.front(), set_indices, ~0u)) {
            add = false;
            break;
          }
        }
        new_output.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
#endif
        ++parent_idx;
      }

      if (add) {
        ctx->outputs->emplace_back(std::move(new_output));
        output_num++;
      }
    }
    return output_num;
  }

  template <QueryType profile>
  uint32_t fromSetNeighborStrategy(QueryVertexID min_parent, uint32_t min_parent_idx, TraverseContext* ctx) const {
    unordered_set<VertexID> visited;
    auto& input = ctx->getCurrentInput();
    uint32_t output_num = 0;
    const auto& parent_match = input.getSet(query_vertex_indices_[min_parent]);
    auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);

    for (VertexID vid : *parent_match) {
      auto out_neighbors = ((G*)(ctx->current_data_graph))->getOutNeighborsWithHint(vid, ALL_LABEL, min_parent_idx);
      for (VertexID key_vertex_id : out_neighbors) {
        if (exceptions.count(key_vertex_id)) {
          continue;
        }
        if (visited.insert(key_vertex_id).second) {
          if (!isInCandidates(key_vertex_id)) {
            continue;
          }
          auto key_out_neighbors =
              ((G*)(ctx->current_data_graph))->getInNeighborsWithHint(key_vertex_id, ALL_LABEL, min_parent_idx);
          // TODO(by) hash key_out_neighbors
          CompressedSubgraphs new_output(input, key_vertex_id, same_label_set_indices_, set_pruning_threshold_);
          if (new_output.empty()) {
            continue;
          }
          bool add = true;
#ifndef USE_FILTER
          // TODO(tatiana): `ExpandSetToKey` requires different groups of same-label indices for the parent sets
          // for active pruning, should use same-label set indices >>>
          unordered_set<uint32_t> set_indices;
          for (uint32_t i = 0; i < new_output.getNumSets(); ++i) {
            set_indices.insert(i);
          }
// <<< for active pruning, should use same-label set indices
#endif
          uint32_t parent_idx = 0;
          for (uint32_t set_vid : parents_) {
            std::vector<VertexID> new_set;
            uint32_t id = query_vertex_indices_[set_vid];
            if
              constexpr(
                  !std::is_same<G, Graph>::value) {  // if using graph view, select neighbors from the right graph part
                key_out_neighbors =
                    ((G*)(ctx->current_data_graph))->getInNeighborsWithHint(key_vertex_id, ALL_LABEL, parent_idx);
              }
            intersect(*input.getSet(id), key_out_neighbors, &new_set);
            if
              constexpr(isProfileMode(profile)) {
                updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.size(), new_set.size());
              }
            if (new_set.empty()) {
              add = false;
              break;
            }
#ifdef USE_FILTER
            new_output.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
            if (filter(new_output)) {
              add = false;
              break;
            }
#else
            if (new_set.size() == 1) {  // actively prune existing sets
              auto pruning_set_indices = set_indices;
              pruning_set_indices.erase(id);
              if (new_output.pruneExistingSets(new_set.front(), pruning_set_indices, ~0u)) {
                add = false;
                break;
              }
            }
            new_output.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
#endif
            ++parent_idx;
          }

          if (add) {
            ctx->outputs->emplace_back(std::move(new_output));
            output_num++;
          }
        }
      }
    }
    return output_num;
  }

  /** Calculate the ideal si count as if we expand parent 1, 2, ... n for n = parents_.size() in normal backtracing
   * implementation. */
  void updateDistinctSICount(ExpandVertexTraverseContext* ctx) const {
    auto& input = ctx->getCurrentInput();
    std::vector<std::vector<VertexID>*> parent_set_ptrs;
    parent_set_ptrs.reserve(parents_.size());
    for (auto parent : parents_) {
      parent_set_ptrs.push_back(input.getSet(query_vertex_indices_[parent]).get());
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
        ctx->updateDistinctSICount(depth, parent_tuple, depth);
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

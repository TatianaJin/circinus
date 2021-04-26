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

#include <string>
#include <tuple>
#include <vector>

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

template <typename G>
class ExpandSetToKeyVertexOperator : public ExpandVertexOperator {
  std::vector<unordered_set<std::string>> parent_tuple_sets_;

 public:
  ExpandSetToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               SubgraphFilter* subgraph_filter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, subgraph_filter) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<QueryType::Execute>(outputs, batch_size);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size,
                                 uint32_t query_type) override {
    if (query_type == 1) {
      return expandInner<QueryType::Profile>(outputs, batch_size);
    }
    CHECK_EQ(query_type, 2) << "unknown query type " << query_type;
    return expandInner<QueryType::ProfileWithMiniIntersection>(outputs, batch_size);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandSetToKeyVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new ExpandSetToKeyVertexOperator(*this);
  }

 protected:
  uint32_t profile() {
    auto& input = (*current_inputs_)[input_index_];
    uint32_t idx = 0;
    for (auto par : parents_) {
      uint32_t loop_num = 0;
      const auto& parent_match = input.getSet(query_vertex_indices_[par]);
      for (VertexID vid : *parent_match) {
        const auto& out_neighbors = ((G*)current_data_graph_)->getOutNeighbors(vid, 0, idx);
        for (uint32_t i = 0; i < out_neighbors.second; ++i) {
          loop_num++;
        }
      }
      DLOG(INFO) << "set out_neighbors num " << loop_num;
      ++idx;
    }
    DLOG(INFO) << "candidates num " << candidates_->size();
    return 0;
  }

  bool isInCandidates(VertexID key) {
    auto lb = std::lower_bound(candidates_->begin(), candidates_->end(), key);
    return lb != candidates_->end() && *lb == key;
  }

  std::tuple<uint32_t, uint32_t, uint32_t> getMinimumParent() {
    uint32_t parent = 0, size = 0xFFFFFFFF, parent_idx = 0;
    auto& input = (*current_inputs_)[input_index_];
    uint32_t idx = 0;
    for (auto par : parents_) {
      auto current_set = input.getSet(query_vertex_indices_[par]);
      uint32_t current_size = 0;
      for (auto set_vertex_id : *current_set) {
        current_size += ((G*)current_data_graph_)->getVertexOutDegree(set_vertex_id, 0, idx);
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
  inline uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    if
      constexpr(isProfileWithMiniIntersectionMode(profile)) { parent_tuple_sets_.resize(parents_.size()); }
    uint32_t output_num = 0;
    while (input_index_ < current_inputs_->size()) {
      auto[min_parent_set_size, min_parent_vertex, min_parent_idx] = getMinimumParent();
      if
        constexpr(isProfileWithMiniIntersectionMode(profile)) { updateDistinctSICount(); }
      if (min_parent_set_size < candidates_->size()) {
        DLOG(INFO) << "fromSetNeighborStrategy";
        output_num += fromSetNeighborStrategy<profile>(outputs, min_parent_vertex, min_parent_idx);
      } else {
        DLOG(INFO) << "fromCandidateStrategy";
        output_num += fromCandidateStrategy<profile>(outputs);
      }
      if
        constexpr(isProfileMode(profile)) {
          total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
        }
      input_index_++;
      if (output_num >= batch_size) {
        break;
      }
    }
    return output_num;
  }

  template <QueryType profile>
  uint32_t fromCandidateStrategy(std::vector<CompressedSubgraphs>* outputs) {
    auto& input = (*current_inputs_)[input_index_];
    uint32_t output_num = 0;
    auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
    for (VertexID key_vertex_id : *candidates_) {
      if (exceptions.count(key_vertex_id)) {
        continue;
      }
      auto key_out_neighbors = ((G*)current_data_graph_)->getOutNeighbors(key_vertex_id, 0, 0);
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
            key_out_neighbors = ((G*)current_data_graph_)->getOutNeighbors(key_vertex_id, 0, parent_idx);
          }
        intersect(*input.getSet(id), key_out_neighbors, &new_set);  // No need for exceptions
        if
          constexpr(isProfileMode(profile)) {
            updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.second, new_set.size());
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
        outputs->emplace_back(std::move(new_output));
        output_num++;
      }
    }
    return output_num;
  }

  template <QueryType profile>
  uint32_t fromSetNeighborStrategy(std::vector<CompressedSubgraphs>* outputs, QueryVertexID min_parent,
                                   uint32_t min_parent_idx) {
    unordered_set<VertexID> visited;
    auto& input = (*current_inputs_)[input_index_];
    uint32_t output_num = 0;
    const auto& parent_match = input.getSet(query_vertex_indices_[min_parent]);
    auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);

    for (VertexID vid : *parent_match) {
      const auto& out_neighbors = ((G*)current_data_graph_)->getOutNeighbors(vid, 0, min_parent_idx);
      for (uint32_t i = 0; i < out_neighbors.second; ++i) {
        VertexID key_vertex_id = out_neighbors.first[i];
        if (exceptions.count(key_vertex_id)) {
          continue;
        }
        if (visited.insert(key_vertex_id).second) {
          if (!isInCandidates(key_vertex_id)) {
            continue;
          }
          auto key_out_neighbors = ((G*)current_data_graph_)->getOutNeighbors(key_vertex_id, 0, min_parent_idx);
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
                key_out_neighbors = ((G*)current_data_graph_)->getOutNeighbors(key_vertex_id, 0, parent_idx);
              }
            intersect(*input.getSet(id), key_out_neighbors, &new_set);
            if
              constexpr(isProfileMode(profile)) {
                updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.second, new_set.size());
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
            outputs->emplace_back(std::move(new_output));
            output_num++;
          }
        }
      }
    }
    return output_num;
  }

  /** Calculate the ideal si count as if we expand parent 1, 2, ... n for n = parents_.size() in normal backtracing
   * implementation. */
  void updateDistinctSICount() {
    auto& input = (*current_inputs_)[input_index_];
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
        distinct_intersection_count_ +=
            parent_tuple_sets_[depth].emplace((char*)parent_tuple.data(), (depth + 1) * sizeof(VertexID)).second;
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

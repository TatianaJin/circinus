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

#include <chrono>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "ops/traverse_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

template <typename G>
class ExpandIntoOperator : public TraverseOperator {
  std::vector<QueryVertexID> parents_;
  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;
  // for profiling
  std::vector<QueryVertexID> key_parents_;  // the key parents of the previous operator
  std::vector<unordered_set<std::string>> parent_tuple_sets_;

 public:
  ExpandIntoOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                     const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                     const std::vector<QueryVertexID>& prev_key_parents, SubgraphFilter* subgraph_filter = nullptr)
      : TraverseOperator(target_vertex, subgraph_filter),
        parents_(parents),
        query_vertex_indices_(query_vertex_indices),
        key_parents_(prev_key_parents) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<QueryType::Execute>(outputs, batch_size);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size,
                                 uint32_t query_type) override {
    if (query_type == 1) {
      return expandInner<QueryType::Profile>(outputs, batch_size);
    } else {
      CHECK_EQ(query_type, 2) << "unknown query type " << query_type;
      return expandInner<QueryType::ProfileWithMiniIntersection>(outputs, batch_size);
    }
  }

  std::vector<BipartiteGraph*> computeBipartiteGraphs(
      const Graph* g, const std::vector<std::vector<VertexID>>& candidate_sets) override {
    std::vector<BipartiteGraph*> ret;
    ret.reserve(parents_.size());
    for (auto parent_vertex : parents_) {
      ret.emplace_back(new BipartiteGraph(parent_vertex, target_vertex_));
      ret.back()->populateGraph(g, &candidate_sets);
    }
    return ret;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandIntoOperator";
    for (auto parent : parents_) {
      DCHECK_EQ(query_vertex_indices_.count(parent), 1);
      ss << ' ' << parent;
    }
    DCHECK_EQ(query_vertex_indices_.count(target_vertex_), 1);
    ss << " -> " << target_vertex_;
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    auto ret = new ExpandIntoOperator(*this);
    return ret;
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    // for profiling distinct si count, consider the cost as expanding from all parents of the current target without
    // compression
    std::vector<VertexID> parent_tuple;
    if
      constexpr(isProfileMode(profile)) {
        parent_tuple.resize(key_parents_.size() + parents_.size());
        parent_tuple_sets_.resize(parents_.size());
      }
    while (input_index_ < current_inputs_->size()) {
      auto input = (*current_inputs_)[input_index_];
      auto key_vertex_id = input.getKeyVal(query_vertex_indices_[target_vertex_]);
      auto key_out_neighbors = ((G*)current_data_graph_)->getOutNeighborsWithHint(key_vertex_id, 0, 0);
      bool add = true;
      if
        constexpr(isProfileMode(profile)) {
          total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
          if
            constexpr(isProfileWithMiniIntersectionMode(profile)) {
              unordered_set<VertexID> prefix_set;
              for (uint32_t i = 0; i < key_parents_.size(); ++i) {
                parent_tuple[i] = input.getKeyVal(query_vertex_indices_[key_parents_[i]]);
                prefix_set.insert(parent_tuple[i]);
              }
              std::vector<std::vector<VertexID>*> parent_set_ptrs;
              parent_set_ptrs.reserve(parents_.size());
              for (auto parent : parents_) {
                parent_set_ptrs.push_back(input.getSet(query_vertex_indices_[parent]).get());
              }
              uint32_t depth = 0, last_depth = parents_.size() - 1;
              std::vector<uint32_t> set_index(parents_.size(), 0);
              while (true) {
                while (set_index[depth] < parent_set_ptrs[depth]->size()) {
                  auto parent_vid = (*parent_set_ptrs[depth])[set_index[depth]];
                  if (prefix_set.count(parent_vid)) {
                    ++set_index[depth];
                    continue;
                  }
                  auto pidx = depth + key_parents_.size();
                  parent_tuple[pidx] = parent_vid;
                  distinct_intersection_count_ +=
                      parent_tuple_sets_[depth]
                          .emplace((char*)parent_tuple.data(), (pidx + 1) * sizeof(VertexID))
                          .second;
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
        }
#ifndef USE_FILTER
      // TODO(tatiana): `ExpandInto` requires different groups of same-label indices for the parent sets
      // for active pruning, should use same-label set indices >>>
      unordered_set<uint32_t> set_indices;
      for (uint32_t i = 0; i < input.getNumSets(); ++i) {
        set_indices.insert(i);
      }  // <<< for active pruning, should use same-label set indices
#endif
      // TODO(tatiana): consider sorting of parents in ascending order of set size, for better pruning?
      uint32_t parent_idx = 0;
      for (QueryVertexID vid : parents_) {
        std::vector<VertexID> new_set;
        uint32_t id = query_vertex_indices_[vid];
        if
          constexpr(!std::is_same<G, Graph>::value) {
            key_out_neighbors = ((G*)current_data_graph_)->getOutNeighborsWithHint(key_vertex_id, 0, parent_idx);
          }
        intersect(*(input.getSet(id)), key_out_neighbors, &new_set);
        if
          constexpr(isProfileMode(profile)) {
            updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.second, new_set.size());
          }
        if (new_set.size() == 0) {
          add = false;
          break;
        }
#ifdef USE_FILTER
        // TODO(tatiana): include same-label keys for checking when sets are updated (ExpandInto and ExpandSettoKey)
        input.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
        if (filter(input)) {
          add = false;
          break;
        }
#else
        if (new_set.size() == 1) {  // actively prune existing sets
          auto pruning_set_indices = set_indices;
          pruning_set_indices.erase(id);
          if (input.pruneExistingSets(new_set.front(), pruning_set_indices, set_pruning_threshold_)) {
            add = false;
            break;
          }
        }
        input.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
#endif
        ++parent_idx;
      }

      if (add) {
        outputs->emplace_back(std::move(input));
        output_num++;
      }
      input_index_++;

      if (output_num >= batch_size) {
        break;
      }
    }
    return output_num;
  }
};

}  // namespace circinus

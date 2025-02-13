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
#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

// TODO(tatiana): handle the case when G is not Graph
template <typename G>
class EnumerateKeyExpandToSetOperator : public ExpandVertexOperator {
  const std::vector<QueryVertexID> keys_to_enumerate_;
  std::vector<uint32_t> existing_key_parent_indices_;
  std::vector<int> cover_table_;
  unordered_map<QueryVertexID, uint32_t> enumerate_key_old_indices_;
  std::vector<std::pair<uint32_t, int>> set_old_to_new_pos_;
  bool need_new_input_ = true;
#ifndef USE_FILTER
  // for set pruning, the indices of the sets with the same label as the target in the output
  unordered_set<uint32_t> set_indices_;
#endif
  std::vector<int> enumerated_key_pruning_indices_;

  /* transient */
  std::vector<uint32_t> enumerate_key_idx_;                     // size = keys_to_enumerate_.size();
  std::vector<std::vector<VertexID>*> enumerate_key_pos_sets_;  // size = keys_to_enumerate_.size();
  std::vector<std::vector<VertexID>> target_sets_;  // now we store and reuse the intermediate intersection results
  unordered_set<VertexID> existing_vertices_;
  uint32_t n_exceptions_ = 0;
  CompressedSubgraphs output_;

  /* for profiling */
  std::vector<unordered_set<std::string>> parent_tuple_sets_;
  std::vector<VertexID> parent_tuple_;

 public:
  EnumerateKeyExpandToSetOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                  const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
                                  const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
                                  const std::vector<QueryVertexID>& keys_to_enumerate,
                                  const std::vector<int>& cover_table,
                                  const std::array<std::vector<uint32_t>, 2>& same_label_indices,
                                  std::vector<int>&& enumerated_key_pruning_indices, SubgraphFilter* filter);

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
    ss << "EnumerateKeyExpandToSetOperator";
    toStringInner(ss);
    ss << " (keys to enumerate";
    for (auto v : keys_to_enumerate_) {
      ss << ' ' << v;
    }
    ss << ')';
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new EnumerateKeyExpandToSetOperator(*this);
  }

 private:
  template <QueryType>
  uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size);

  template <QueryType>
  bool expandInner();
};

template <typename G>
EnumerateKeyExpandToSetOperator<G>::EnumerateKeyExpandToSetOperator(
    const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
    const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
    const std::vector<QueryVertexID>& keys_to_enumerate, const std::vector<int>& cover_table,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices, std::vector<int>&& enumerated_key_pruning_indices,
    SubgraphFilter* subgraph_filter)
    : ExpandVertexOperator(parents, target_vertex, output_query_vertex_indices, same_label_indices[1],
                           same_label_indices[0], ~0u, subgraph_filter),
      keys_to_enumerate_(keys_to_enumerate),
      cover_table_(cover_table),
      enumerated_key_pruning_indices_(std::move(enumerated_key_pruning_indices)),
      enumerate_key_idx_(keys_to_enumerate.size(), 0),
      enumerate_key_pos_sets_(keys_to_enumerate.size()),
      target_sets_(keys_to_enumerate.size() + 1),
      output_(output_query_vertex_indices.size() - 1 - output_query_vertex_indices.at(target_vertex),
              output_query_vertex_indices.size()) {
  CHECK(subgraph_filter != nullptr);
  CHECK_EQ(enumerated_key_pruning_indices_.size(), keys_to_enumerate_.size());

  unordered_set<QueryVertexID> keys_to_enumerate_set(keys_to_enumerate.begin(), keys_to_enumerate.end());
  // get existing key parent indices
  existing_key_parent_indices_.reserve(parents_.size() - keys_to_enumerate_.size());
  for (auto v : parents) {
    if (keys_to_enumerate_set.count(v) == 0) {
      existing_key_parent_indices_.push_back(query_vertex_indices_[v]);
    }
  }

#ifndef USE_FILTER
  // get target pruning set indices
  if (same_label_set_indices_.size() - keys_to_enumerate_.size() > 0) {
    auto pruning_sets = subgraph_filter_->getPruningSets(0);
    set_pruning_threshold_ = subgraph_filter_->getSetPruningThreshold(0);
    CHECK(pruning_sets != nullptr);
    set_indices_.insert(pruning_sets->begin(), pruning_sets->end());
  }
#endif
  for (auto v : keys_to_enumerate_) {
    auto pos =
        std::find(same_label_set_indices_.begin(), same_label_set_indices_.end(), input_query_vertex_indices.at(v));
    if (pos != same_label_set_indices_.end()) same_label_set_indices_.erase(pos);
  }

  // get index mapping of enumerated key vertex indices and set vertex indices
  uint32_t n_input_keys = 0;
  for (auto& pair : input_query_vertex_indices) {
    if (cover_table_[pair.first] != 1) {
      auto new_pos = query_vertex_indices_.at(pair.first);
      set_old_to_new_pos_.emplace_back(pair.second, new_pos);
    } else {
      n_input_keys += (keys_to_enumerate_set.count(pair.first) == 0);
    }
  }
  for (auto v : keys_to_enumerate_) {
    DCHECK(input_query_vertex_indices.count(v));
    enumerate_key_old_indices_[v] = input_query_vertex_indices.at(v);
    CHECK_EQ(query_vertex_indices_[v], n_input_keys) << v;  // assume contiguous indices after existing keys
    ++n_input_keys;
  }
  // assume target is the last set in output
  CHECK_EQ(query_vertex_indices_[target_vertex_], output_query_vertex_indices.size() - n_input_keys - 1);
}

template <typename G>
template <QueryType profile>
uint32_t EnumerateKeyExpandToSetOperator<G>::expandInner(std::vector<CompressedSubgraphs>* outputs,
                                                         uint32_t batch_size) {
  uint32_t n_outputs = 0;
  if
    constexpr(isProfileWithMiniIntersectionMode(profile)) {
      parent_tuple_sets_.resize(parents_.size());
      parent_tuple_.resize(parents_.size());
    }
  while (input_index_ < current_inputs_->size()) {
    if (need_new_input_) {
      // find next input with non-empty candidate target set
      while (input_index_ < current_inputs_->size() && expandInner<profile>()) {
        total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
        ++input_index_;
      }
      if (input_index_ == current_inputs_->size()) {  // all inputs are consumed
        return n_outputs;
      }
      total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
      // reset index
      enumerate_key_idx_[0] = 0;
      // get sets to enumerate as compression key
      auto& input = (*current_inputs_)[input_index_];
      uint32_t depth = 0;
      for (auto v : keys_to_enumerate_) {
        auto pos = enumerate_key_old_indices_.at(v);
        enumerate_key_pos_sets_[depth] = input.getSet(pos).get();
        ++depth;
      }
      // set existing matched vertices in output
      std::copy(input.getKeys().begin(), input.getKeys().end(), output_.getKeys().begin());
      DCHECK_EQ(set_old_to_new_pos_.size(), output_.getNumSets() - 1);
      for (auto& pair : set_old_to_new_pos_) {
        output_.UpdateSets(pair.second, input.getSet(pair.first));
      }
      need_new_input_ = false;
    }
    const auto& input = (*current_inputs_)[input_index_];
    const uint32_t enumerate_key_size = keys_to_enumerate_.size();
    uint32_t enumerate_key_depth = existing_vertices_.size() - n_exceptions_;
    if (enumerate_key_depth != 0) {  // if continuing from a previously processed input
      DCHECK_EQ(enumerate_key_depth, enumerate_key_size - 1)
          << "existing_vertices_.size()=" << existing_vertices_.size() << ", input.getNumKeys()=" << input.getNumKeys()
          << '/' << n_exceptions_;
    }

    while (true) {
      while (enumerate_key_idx_[enumerate_key_depth] < enumerate_key_pos_sets_[enumerate_key_depth]->size()) {
        // update target set
        auto key_vid = (*enumerate_key_pos_sets_[enumerate_key_depth])[enumerate_key_idx_[enumerate_key_depth]];
        if (existing_vertices_.count(key_vid) == 1) {  // skip current key to ensure isomorphism
          ++enumerate_key_idx_[enumerate_key_depth];
          continue;
        }
        DCHECK(target_sets_[enumerate_key_depth + 1].empty()) << enumerate_key_depth;
        if (target_sets_[enumerate_key_depth].empty()) {
          DCHECK_EQ(enumerate_key_depth, 0);
          intersect(*candidates_, ((G*)current_data_graph_)->getOutNeighborsWithHint(key_vid, 0, 0),
                    &target_sets_[enumerate_key_depth + 1], existing_vertices_);
        } else {
          intersect(target_sets_[enumerate_key_depth],
                    ((G*)current_data_graph_)->getOutNeighborsWithHint(key_vid, 0, 0),
                    &target_sets_[enumerate_key_depth + 1], existing_vertices_);
        }

        if
          constexpr(isProfileMode(profile)) {
            updateIntersectInfo((target_sets_[enumerate_key_depth].empty() ? candidates_->size()
                                                                           : target_sets_[enumerate_key_depth].size()) +
                                    ((G*)current_data_graph_)->getVertexOutDegreeWithHint(key_vid, 0, 0),
                                target_sets_[enumerate_key_depth + 1].size());
            if
              constexpr(isProfileWithMiniIntersectionMode(profile)) {
                auto pidx = enumerate_key_depth + existing_key_parent_indices_.size();  // parent index
                parent_tuple_[pidx] = key_vid;
                distinct_intersection_count_ +=
                    parent_tuple_sets_[pidx].emplace((char*)parent_tuple_.data(), (pidx + 1) * sizeof(VertexID)).second;
              }
          }
        if (target_sets_[enumerate_key_depth + 1].empty()) {
          ++enumerate_key_idx_[enumerate_key_depth];
          continue;
        }
        if (enumerate_key_depth == enumerate_key_size - 1) {  // the last key query vertex to enumerate, ready to output
          auto& target_set = target_sets_.back();
          auto output = output_;
          // set the enumerated keys in the output
          bool skip = false;
          for (uint32_t key_i = 0; key_i < enumerate_key_size; ++key_i) {
            auto key = (*enumerate_key_pos_sets_[key_i])[enumerate_key_idx_[key_i]];
            output.UpdateKey(input.getNumKeys() + key_i, key);
            if (enumerated_key_pruning_indices_[key_i] != -1) {
              const auto& pruning_sets = *subgraph_filter_->getPruningSets(enumerated_key_pruning_indices_[key_i]);
              unordered_set<uint32_t> indices(pruning_sets.begin(), pruning_sets.end());
              indices.erase(output.getNumSets() - 1);
              auto thres = subgraph_filter_->getSetPruningThreshold(enumerated_key_pruning_indices_[key_i]);
#ifdef USE_FILTER
              bool recursive_prune = false;  // only prune by key
#else
              bool recursive_prune = true;  // recursively prune
#endif
              if (output.pruneExistingSets(key, indices, thres, recursive_prune)) {
                skip = true;
                break;
              }
            }
          }
          if (skip) {
            target_set.clear();
            ++enumerate_key_idx_[enumerate_key_depth];
            continue;
          }
#ifdef USE_FILTER
          // set the target set
          output.UpdateSets(output.getNumSets() - 1, std::make_shared<std::vector<VertexID>>(std::move(target_set)));
          if (filter(output)) {
            ++enumerate_key_idx_[enumerate_key_depth];
            continue;
          }
#else
          // set the target set
          if (target_set.size() == 1) {
            auto indices = set_indices_;
            if (output.pruneExistingSets(target_set.front(), indices, set_pruning_threshold_)) {
              ++enumerate_key_idx_[enumerate_key_depth];
              continue;
            }
          }
          output.UpdateSets(output.getNumSets() - 1, std::make_shared<std::vector<VertexID>>(std::move(target_set)));
#endif

          outputs->push_back(std::move(output));
          ++enumerate_key_idx_[enumerate_key_depth];
          if (++n_outputs == batch_size) {
            return n_outputs;
          }
        } else {  // dfs next key depth
          existing_vertices_.insert(key_vid);
          ++enumerate_key_depth;
          enumerate_key_idx_[enumerate_key_depth] = 0;  // start from the first match in the next depth
        }
      }
      target_sets_[enumerate_key_depth].clear();
      if (enumerate_key_depth == 0) {  // current input is completely processed
        need_new_input_ = true;
        ++input_index_;
        break;
      }
      --enumerate_key_depth;
      existing_vertices_.erase(
          (*enumerate_key_pos_sets_[enumerate_key_depth])[enumerate_key_idx_[enumerate_key_depth]]);
      ++enumerate_key_idx_[enumerate_key_depth];
    }
  }
  return n_outputs;
}

template <typename G>
template <QueryType profile>
bool EnumerateKeyExpandToSetOperator<G>::expandInner() {  // handles a new input and init the transient states
  auto& input = (*current_inputs_)[input_index_];
  auto& target_set = target_sets_.front();
  existing_vertices_ = input.getKeyMap();
  input.getExceptions(existing_vertices_, {}, same_label_set_indices_);
  n_exceptions_ = existing_vertices_.size();

  for (uint32_t i = 0; i < existing_key_parent_indices_.size(); ++i) {
    uint32_t key_vid = input.getKeyVal(existing_key_parent_indices_[i]);
    auto target_set_size = target_set.size();
    (void)target_set_size;
    if (i == 0) {
      intersect(*candidates_, ((G*)current_data_graph_)->getOutNeighborsWithHint(key_vid, 0, 0), &target_set,
                existing_vertices_);
    } else {
      intersectInplace(target_set, ((G*)current_data_graph_)->getOutNeighborsWithHint(key_vid, 0, 0), &target_set);
    }
    if
      constexpr(isProfileMode(profile)) {
        updateIntersectInfo(target_set_size + ((G*)current_data_graph_)->getVertexOutDegreeWithHint(key_vid, 0, 0),
                            target_set.size());
        if
          constexpr(isProfileWithMiniIntersectionMode(profile)) {
            parent_tuple_[i] = key_vid;
            distinct_intersection_count_ +=
                parent_tuple_sets_[i].emplace((char*)parent_tuple_.data(), (i + 1) * sizeof(VertexID)).second;
          }
      }
    if (target_set.empty()) {
      return true;
    }
  }

  return candidates_->empty() || (target_set.empty() && existing_key_parent_indices_.size() > 0);
}

}  // namespace circinus

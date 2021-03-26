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

#include "ops/enumerate_key_expand_to_set_operator.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

EnumerateKeyExpandToSetOperator::EnumerateKeyExpandToSetOperator(
    const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
    const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
    const std::vector<QueryVertexID>& keys_to_enumerate, const std::vector<int>& cover_table)
    : ExpandVertexOperator(parents, target_vertex, output_query_vertex_indices),
      keys_to_enumerate_(keys_to_enumerate),
      cover_table_(cover_table),
      enumerate_key_idx_(keys_to_enumerate.size(), 0),
      enumerate_key_pos_sets_(keys_to_enumerate.size()),
      target_sets_(keys_to_enumerate.size() + 1) {
  existing_key_parents_.reserve(parents_.size() - keys_to_enumerate_.size());
  unordered_set<QueryVertexID> keys_to_enumerate_set(keys_to_enumerate.begin(), keys_to_enumerate.end());
  for (auto v : parents) {
    if (keys_to_enumerate_set.count(v) == 0) {
      existing_key_parents_.push_back(v);
    }
  }
  uint32_t n_input_keys = 0;
  for (auto& pair : input_query_vertex_indices) {
    if (cover_table_[pair.first] != 1) {
      set_old_to_new_pos_.emplace_back(pair.second, query_vertex_indices_.at(pair.first));
    } else {
      n_input_keys += (keys_to_enumerate_set.count(pair.first) == 0);
    }
  }
  n_input_keys_ = n_input_keys;
  for (auto v : keys_to_enumerate_) {
    DCHECK(input_query_vertex_indices.count(v));
    enumerate_key_old_indices_[v] = input_query_vertex_indices.at(v);
    CHECK_EQ(query_vertex_indices_[v], n_input_keys);  // assume contiguous indices after existing keys
    ++n_input_keys;
  }
  // assume target is the last set in output
  CHECK_EQ(query_vertex_indices_[target_vertex_], output_query_vertex_indices.size() - n_input_keys - 1);
}

template <QueryType profile>
uint32_t EnumerateKeyExpandToSetOperator::expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
  uint32_t n_outputs = 0;
  if
    constexpr(isProfileMode(profile)) {
      parent_tuple_sets_.resize(parents_.size());
      parent_tuple_.resize(parents_.size());
    }

  while (input_index_ < current_inputs_->size()) {
    if (target_sets_.front().empty()) {
      // find next input with non-empty candidate target set
      while (input_index_ < current_inputs_->size() && expandInner<profile>()) {
        ++input_index_;
      }
      if (target_sets_.front().empty()) {  // all inputs are consumed
        return 0;
      }
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
      existing_key_vertices_.clear();
      existing_key_vertices_.insert(input.getKeys().begin(), input.getKeys().end());
      DCHECK_EQ(existing_key_vertices_.size(), n_input_keys_) << input.getNumKeys();
    }

    const auto& input = (*current_inputs_)[input_index_];
    const uint32_t enumerate_key_size = keys_to_enumerate_.size();
    // set existing matched vertices in output
    CompressedSubgraphs output(n_input_keys_ + enumerate_key_size, input.getNumVertices() + 1);
    std::copy(input.getKeys().begin(), input.getKeys().end(), output.getKeys().begin());
    DCHECK_EQ(set_old_to_new_pos_.size(), output.getNumSets() - 1);
    for (auto& pair : set_old_to_new_pos_) {
      output.UpdateSets(pair.second, input.getSet(pair.first));
    }  // TODO(tatiana): partially set output
    uint32_t enumerate_key_depth = existing_key_vertices_.size() - n_input_keys_;
    if (enumerate_key_depth != 0) {  // if continuing from a previously processed input
      CHECK_EQ(enumerate_key_depth, enumerate_key_size - 1)
          << "existing_key_vertices_.size()=" << existing_key_vertices_.size()
          << ", input.getNumKeys()=" << input.getNumKeys() << '/' << n_input_keys_;
    }
    while (true) {
      while (enumerate_key_idx_[enumerate_key_depth] < enumerate_key_pos_sets_[enumerate_key_depth]->size()) {
        // update target set
        auto key_vid = (*enumerate_key_pos_sets_[enumerate_key_depth])[enumerate_key_idx_[enumerate_key_depth]];
        if (existing_key_vertices_.count(key_vid) == 1) {  // skip current key to ensure isomorphism
          ++enumerate_key_idx_[enumerate_key_depth];
          continue;
        }
        DCHECK(target_sets_[enumerate_key_depth + 1].empty());
        intersect(target_sets_[enumerate_key_depth], current_data_graph_->getOutNeighbors(key_vid),
                  &target_sets_[enumerate_key_depth + 1]);
        if
          constexpr(isProfileMode(profile)) {
            updateIntersectInfo(
                target_sets_[enumerate_key_depth].size() + current_data_graph_->getVertexOutDegree(key_vid),
                target_sets_[enumerate_key_depth + 1].size());
            auto pidx = enumerate_key_depth + existing_key_parents_.size();  // parent index
            parent_tuple_[pidx] = key_vid;
            distinct_intersection_count_ +=
                parent_tuple_sets_[pidx].emplace((char*)parent_tuple_.data(), (pidx + 1) * sizeof(VertexID)).second;
          }
        if (target_sets_[enumerate_key_depth + 1].empty()) {
          ++enumerate_key_idx_[enumerate_key_depth];
          continue;
        }
        if (enumerate_key_depth == enumerate_key_size - 1) {
          // the last key query vertex to enumerate, ready to output
          auto& target_set = target_sets_.back();
          if (!target_set.empty()) {
            for (uint32_t key_i = 0; key_i < enumerate_key_size; ++key_i) {
              output.UpdateKey(input.getNumKeys() + key_i,
                               (*enumerate_key_pos_sets_[key_i])[enumerate_key_idx_[key_i]]);
            }
            output.UpdateSets(output.getNumSets() - 1, std::make_shared<std::vector<VertexID>>(std::move(target_set)));
            outputs->push_back(output);
            ++enumerate_key_idx_[enumerate_key_depth];
            if (++n_outputs == batch_size) {
              return n_outputs;
            }
          } else {
            ++enumerate_key_idx_[enumerate_key_depth];
          }
        } else {  // dfs next key depth
          existing_key_vertices_.insert(key_vid);
          ++enumerate_key_depth;
          enumerate_key_idx_[enumerate_key_depth] = 0;  // start from the first match in the next depth
        }
      }
      target_sets_[enumerate_key_depth].clear();
      if (enumerate_key_depth == 0) {  // current input is completely processed
        ++input_index_;
        break;
      }
      --enumerate_key_depth;
      existing_key_vertices_.erase(
          (*enumerate_key_pos_sets_[enumerate_key_depth])[enumerate_key_idx_[enumerate_key_depth]]);
      ++enumerate_key_idx_[enumerate_key_depth];
    }
  }
  return n_outputs;
}

std::string EnumerateKeyExpandToSetOperator::toString() const {
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

template <QueryType profile>
bool EnumerateKeyExpandToSetOperator::expandInner() {  // handles a new input and init the transient states
  auto& input = (*current_inputs_)[input_index_];
  auto& target_set = target_sets_.front();
  target_set = *candidates_;
  for (uint32_t i = 0; i < existing_key_parents_.size(); ++i) {
    uint32_t key = query_vertex_indices_[existing_key_parents_[i]];
    uint32_t key_vid = input.getKeyVal(key);
    auto target_set_size = target_set.size();
    intersectInplace(target_set, current_data_graph_->getOutNeighbors(key_vid), &target_set);
    if
      constexpr(isProfileMode(profile)) {
        updateIntersectInfo(target_set_size + current_data_graph_->getVertexOutDegree(key_vid), target_set.size());
        parent_tuple_[i] = key_vid;
        distinct_intersection_count_ +=
            parent_tuple_sets_[i].emplace((char*)parent_tuple_.data(), (i + 1) * sizeof(VertexID)).second;
      }

    if (target_set.empty()) {
      break;
    }
  }
  target_set.erase(std::remove_if(target_set.begin(), target_set.end(),
                                  [&input](VertexID set_vertex) { return input.isExisting(set_vertex); }),
                   target_set.end());
  return target_set.empty();
}

}  // namespace circinus

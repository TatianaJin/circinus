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
#include <cinttypes>
#include <memory>
#include <string>
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

// linear scan based implementation
template <class ForwardIt, typename T>
inline ForwardIt lower_bound(ForwardIt first, ForwardIt last, const T& value) {
  while (first != last && *first < value) {
    ++first;
  }
  return first;
}

/** A group of subgraphs compressed by some key vertices, which constitute a vertex cover for all subgraphs. */
class CompressedSubgraphs {
  std::vector<VertexID> keys_;
  std::vector<VertexSet> sets_;

 public:
  /**
   * @param key_size The number of key vertices.
   * @param n_vertices The number of vertices in each compressed subgraph.
   */
  CompressedSubgraphs(uint32_t key_size, uint32_t n_vertices) : keys_(key_size), sets_(n_vertices - key_size) {}

  /** Construct a CompressedSubgraphs that contains only a single-vertex subgraph.
   * @param key The vertex id.
   */
  explicit CompressedSubgraphs(VertexID key) : keys_(1) { keys_.front() = key; }

  /** Construct a CompressedSubgraphs that contains only a single-vertex subgraph.
   * @param key The vertex id.
   */
  explicit CompressedSubgraphs(VertexSet&& set) : sets_(1) { sets_.front() = std::move(set); }

  /**
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param key The new vertex expanded in this CompressedSubgraphs, which is a key.
   * @param pruning_set_indices The indices of the existing sets to prune.
   * @param set_pruning_threshold The existing sets are only pruned if their sizes are smaller than the threshold.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, VertexID key,
                      const std::vector<uint32_t>& pruning_set_indices, uint64_t set_pruning_threshold,
                      bool recursive_prune = true)
      : keys_(subgraphs.getNumKeys() + 1), sets_(subgraphs.getNumSets()) {
    unordered_set<uint32_t> set_indices(pruning_set_indices.begin(), pruning_set_indices.end());
    sets_ = subgraphs.sets_;
    // isomorphism check for pruning sets that are completely conflicting with keys
    if (pruneExistingSets(key, set_indices, set_pruning_threshold, recursive_prune)) {
      keys_.clear();
      return;
    }
    std::copy(subgraphs.keys_.begin(), subgraphs.keys_.end(), keys_.begin());
    keys_.back() = key;
  }

  /**
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param new_set The set of new vertices expanded in this CompressedSubgraphs, which is not in key.
   * @param pruning_set_indices The indices of the existing sets to prune.
   * @param set_pruning_threshold The existing sets are only pruned if their sizes are smaller than the threshold.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, std::vector<VertexID>&& new_set,
                      const std::vector<uint32_t>& pruning_set_indices, uint64_t set_pruning_threshold)
      : keys_(subgraphs.getNumKeys()), sets_(subgraphs.getNumSets() + 1) {
    std::copy(subgraphs.sets_.begin(), subgraphs.sets_.end(), sets_.begin());
    if (new_set.size() == 1) {
      unordered_set<uint32_t> set_indices(pruning_set_indices.begin(), pruning_set_indices.end());
      if (pruneExistingSets(new_set.front(), set_indices, set_pruning_threshold)) {
        keys_.clear();
        sets_.clear();
        return;
      }
    }
    keys_ = subgraphs.keys_;
    sets_.back() = std::make_shared<std::vector<VertexID>>(std::move(new_set));
  }

  /** This constructor is used when expanding from a non-key vertex to a new key vertex, and the matches of the parent
   * query vertex is to be regrouped due to adding the new key.
   *
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param replacing_set_index The index of the set to be replaced in `subgraphs`.
   * @param new_set The set of new vertices to be put at `replacing_set_index`.
   * @param key The new vertex expanded in this CompressedSubgraphs, which is a key.
   * @param pruning_set_indices The indices of the existing sets to prune.
   * @param set_pruning_threshold The existing sets are only pruned if their sizes are smaller than the threshold.
   * @param prune_by_set Whether to actively prune when the set size is 1.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, uint32_t replacing_set_index, VertexSet&& new_set,
                      VertexID key, const std::vector<uint32_t>& pruning_set_indices, uint64_t set_pruning_threshold,
                      bool prune_by_set)
      : keys_(subgraphs.getNumKeys() + 1), sets_(subgraphs.getNumSets()) {
    sets_ = subgraphs.sets_;
    unordered_set<uint32_t> set_indices(pruning_set_indices.begin(), pruning_set_indices.end());
    if (pruneExistingSets(key, set_indices, set_pruning_threshold)) {
      // TODO(tatiana): set label may not be the same as the new key
      // (prune_by_set && new_set->size() == 1 && pruneExistingSets(new_set->front(), set_indices,
      // set_pruning_threshold))
      keys_.clear();
      sets_.clear();
      return;
    }
    std::copy(subgraphs.keys_.begin(), subgraphs.keys_.end(), keys_.begin());
    keys_.back() = key;
    sets_[replacing_set_index] = std::move(new_set);
  }

  /**
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param key The new vertex expanded in this CompressedSubgraphs, which is a key.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, VertexID key)
      : keys_(subgraphs.getNumKeys() + 1), sets_(subgraphs.getNumSets()) {
    sets_ = subgraphs.sets_;
    std::copy(subgraphs.keys_.begin(), subgraphs.keys_.end(), keys_.begin());
    keys_.back() = key;
  }

  /**
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param new_set The set of new vertices expanded in this CompressedSubgraphs, which is not in key.
   * @param pruning_set_indices The indices of the existing sets to prune.
   * @param set_pruning_threshold The existing sets are only pruned if their sizes are smaller than the threshold.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, std::vector<VertexID>&& new_set)
      : keys_(subgraphs.getNumKeys()), sets_(subgraphs.getNumSets() + 1) {
    keys_ = subgraphs.keys_;
    std::copy(subgraphs.sets_.begin(), subgraphs.sets_.end(), sets_.begin());
    sets_.back() = std::make_shared<std::vector<VertexID>>(std::move(new_set));
  }

  /** This constructor is used when expanding from a non-key vertex to a new key vertex, and the matches of the parent
   * query vertex is to be regrouped due to adding the new key.
   *
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param replacing_set_index The index of the set to be replaced in `subgraphs`.
   * @param new_set The set of new vertices to be put at `replacing_set_index`.
   * @param key The new vertex expanded in this CompressedSubgraphs, which is a key.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, uint32_t replacing_set_index, VertexSet&& new_set,
                      VertexID key)
      : keys_(subgraphs.getNumKeys() + 1), sets_(subgraphs.getNumSets()) {
    std::copy(subgraphs.keys_.begin(), subgraphs.keys_.end(), keys_.begin());
    keys_.back() = key;
    sets_ = subgraphs.sets_;
    sets_[replacing_set_index] = std::move(new_set);
  }

  inline size_t getNumKeys() const { return keys_.size(); }
  inline size_t getNumSets() const { return sets_.size(); }
  inline size_t getNumVertices() const { return getNumKeys() + getNumSets(); }
  inline uint64_t getNumSubgraphs() const {
    if (sets_.empty()) {
      return !keys_.empty();
    }
    uint64_t n_subgraphs = sets_.front()->size();
    for (uint32_t i = 1; i < sets_.size(); ++i) {
      n_subgraphs *= sets_[i]->size();
    }
    return n_subgraphs;
  }

  uint64_t getNumIsomorphicSubgraphs(uint64_t limit = ~0u) const {
    if (sets_.empty()) {
      return !keys_.empty();
    }
    auto sets_sorted_by_size = sets_;
    std::sort(sets_sorted_by_size.begin(), sets_sorted_by_size.end(),
              [](const auto& set1, const auto& set2) { return set1->size() < set2->size(); });
    // dfs sets_ chain
    uint64_t count = 0;
    std::vector<uint32_t> set_index(sets_.size(), 0);
    unordered_set<VertexID> existing_vertices;
    existing_vertices.reserve(getNumVertices());
    existing_vertices.insert(keys_.begin(), keys_.end());
    uint32_t last_depth = sets_.size() - 1;
    uint32_t current_depth = 0;
    while (true) {
      while (set_index[current_depth] < (*sets_sorted_by_size[current_depth]).size()) {
        auto v = (*sets_sorted_by_size[current_depth])[set_index[current_depth]];
        ++set_index[current_depth];
        if (existing_vertices.count(v) == 0) {  // v is valid
          if (current_depth == last_depth) {    // reaching a leave in dfs
            if (++count == limit) return count;
          } else {
            existing_vertices.insert(v);
            ++current_depth;
            set_index[current_depth] = 0;  // start from the first vertex in the next set
          }
        }
      }
      if (current_depth == 0) {
        break;
      }
      --current_depth;
      existing_vertices.erase((*sets_sorted_by_size[current_depth])[set_index[current_depth] - 1]);
    }
    return count;
  }

  bool isExisting(uint32_t key) const {
    for (uint32_t existing_key : keys_) {
      if (existing_key == key) {
        return true;
      }
    }
    return false;
  }

  std::string toString() const {
    std::string s = "{";
    s += "key_size:" + std::to_string(keys_.size()) + "," + "set_size:" + std::to_string(sets_.size()) + ",";
    for (VertexID key : keys_) {
      s += std::to_string(key) + ",";
    }
    for (auto& set : sets_) {
      s += "[";
      for (auto vid : *set) {
        s += std::to_string(vid) + ",";
      }
      s += "],";
    }
    s += "}";
    return s;
  }

  /** Get the value of the key vertex at key_idx. */
  VertexID getKeyVal(uint32_t key_idx) const { return keys_[key_idx]; }

  const std::vector<VertexID>& getKeys() const { return keys_; }
  std::vector<VertexID>& getKeys() { return keys_; }

  unordered_set<VertexID> getKeyMap() const { return unordered_set<VertexID>(keys_.begin(), keys_.end()); }

  // exceptions include keys and single-element sets with the same label
  void getExceptions(unordered_set<VertexID>& exception, const std::vector<uint32_t>& exception_key_indices,
                     const std::vector<uint32_t>& exception_set_indices) const {
    exception.reserve(exception_key_indices.size() + exception_set_indices.size());
    for (auto idx : exception_key_indices) {
      exception.insert(keys_[idx]);
    }

    for (auto idx : exception_set_indices) {
      auto& set = *sets_[idx];
      if (set.size() == 1) {
        exception.insert(set.front());
      }
    }
  }

  // exceptions include keys and single-element sets with the same label
  unordered_set<VertexID> getExceptions(const std::vector<uint32_t>& exception_key_indices,
                                        const std::vector<uint32_t>& exception_set_indices) const {
    unordered_set<VertexID> exception;
    getExceptions(exception, exception_key_indices, exception_set_indices);
    return exception;
  }

  /** Get the matching set of the non-key vertex at key_idx. */
  const VertexSet& getSet(uint32_t key_idx) const { return sets_[key_idx]; }

  void UpdateSets(uint32_t set_idx, VertexSet&& new_set) { sets_[set_idx] = std::move(new_set); }
  void UpdateSets(uint32_t set_idx, const VertexSet& new_set) { sets_[set_idx] = new_set; }

  /** Update the key vertex at key_idx to val. */
  void UpdateKey(uint32_t key_idx, VertexID val) { keys_[key_idx] = val; }
  /** Add a vertex val to the vertex set at set_idx. */
  void UpdateSet(uint32_t set_idx, VertexID val) { sets_[set_idx]->push_back(val); }

  void logString(std::ostream& ss) const {
    for (auto key : keys_) {
      ss << key << ',';
    }
    for (auto& set : sets_) {
      for (auto v : *set) {
        ss << v << ' ';
      }
      ss << ',';
    }
    ss << std::endl;
  }

  bool empty() const { return keys_.empty(); }

  bool pruneExistingSets(VertexID v, unordered_set<uint32_t>& set_indices, uint32_t set_size_threshold,
                         bool recursive_prune = true) {
    unordered_map<VertexID, uint32_t> new_v;  // vertex, set_index
    for (uint32_t i : set_indices) {
      auto& set = *sets_[i];
      if (set.size() <= set_size_threshold) {
        auto lb = circinus::lower_bound(set.begin(), set.end(), v);
        if (lb != set.end() && *lb == v) {  // conflict found
          if (set.size() == 1) {
            return true;
          }
          auto new_set = std::make_shared<std::vector<VertexID>>();
          new_set->insert(new_set->end(), set.begin(), lb);
          new_set->insert(new_set->end(), lb + 1, set.end());
          // recursively prune when set size becomes 1
          if (new_set->size() == 1 && !new_v.insert({new_set->front(), i}).second) {
            return true;
          }
          sets_[i] = std::move(new_set);
        }
      }
    }
    if (new_v.empty()) {
      return false;
    }
    if (recursive_prune) {
      for (auto& pair : new_v) {
        set_indices.erase(pair.second);
        if (pruneExistingSets(pair.first, set_indices, set_size_threshold)) {
          return true;
        }
      }
    }
    return false;
  }
};

}  // namespace circinus

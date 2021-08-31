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
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
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
   * {{key indices},{set indices}}, set indices must not be empty.
   */
  using PruningIndexGroups = std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>;

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
    if (reset(subgraphs, key, pruning_set_indices, set_pruning_threshold, recursive_prune) == nullptr) {
      keys_.clear();
    }
  }

  inline CompressedSubgraphs* reset(const CompressedSubgraphs& subgraphs, VertexID key,
                                    const std::vector<uint32_t>& pruning_set_indices, uint64_t set_pruning_threshold,
                                    bool recursive_prune = true) {
    DCHECK_EQ(sets_.size(), subgraphs.sets_.size());
    DCHECK_EQ(keys_.size(), subgraphs.keys_.size() + 1);
    sets_ = subgraphs.sets_;
    unordered_set<uint32_t> set_indices(pruning_set_indices.begin(), pruning_set_indices.end());
    // isomorphism check for pruning sets that are completely conflicting with keys
    if (pruneExistingSets(key, set_indices, set_pruning_threshold, recursive_prune)) {
      return nullptr;
    }
    std::copy(subgraphs.keys_.begin(), subgraphs.keys_.end(), keys_.begin());
    keys_.back() = key;
    return this;
  }

  /**
   * For pruning using size-one set
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param new_set The set of new vertices expanded in this CompressedSubgraphs, which is not in key.
   * @param pruning_set_indices The indices of the existing sets to prune.
   * @param set_pruning_threshold The existing sets are only pruned if their sizes are smaller than the threshold.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, std::vector<VertexID>&& new_set,
                      const std::vector<uint32_t>& pruning_set_indices, uint64_t set_pruning_threshold)
      : keys_(subgraphs.getNumKeys()), sets_(subgraphs.getNumSets() + 1) {
    if (reset(subgraphs, std::move(new_set), pruning_set_indices, set_pruning_threshold) == nullptr) {
      keys_.clear();
      sets_.clear();
    }
  }

  inline CompressedSubgraphs* reset(const CompressedSubgraphs& subgraphs, std::vector<VertexID>&& new_set,
                                    const std::vector<uint32_t>& pruning_set_indices, uint64_t set_pruning_threshold) {
    std::copy(subgraphs.sets_.begin(), subgraphs.sets_.end(), sets_.begin());
    if (new_set.size() == 1) {
      unordered_set<uint32_t> set_indices(pruning_set_indices.begin(), pruning_set_indices.end());
      if (pruneExistingSets(new_set.front(), set_indices, set_pruning_threshold)) {
        return nullptr;
      }
    }
    keys_ = subgraphs.keys_;
    sets_.back() = std::make_shared<std::vector<VertexID>>(std::move(new_set));
    return this;
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
    if (reset(subgraphs, replacing_set_index, std::move(new_set), key, pruning_set_indices, set_pruning_threshold,
              prune_by_set) == nullptr) {
      keys_.clear();
      sets_.clear();
    }
  }

  inline CompressedSubgraphs* reset(const CompressedSubgraphs& subgraphs, uint32_t replacing_set_index,
                                    VertexSet&& new_set, VertexID key, const std::vector<uint32_t>& pruning_set_indices,
                                    uint64_t set_pruning_threshold, bool prune_by_set) {
    DCHECK_EQ(subgraphs.getNumKeys() + 1, keys_.size());
    sets_ = subgraphs.sets_;
    unordered_set<uint32_t> set_indices(pruning_set_indices.begin(), pruning_set_indices.end());
    if (pruneExistingSets(key, set_indices, set_pruning_threshold)) {
      // TODO(tatiana): set label may not be the same as the new key
      // (prune_by_set && new_set->size() == 1 && pruneExistingSets(new_set->front(), set_indices,
      // set_pruning_threshold))
      return nullptr;
    }
    std::copy(subgraphs.keys_.begin(), subgraphs.keys_.end(), keys_.begin());
    keys_.back() = key;
    sets_[replacing_set_index] = std::move(new_set);
    return this;
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
    reset(subgraphs, std::move(new_set));
  }

  inline CompressedSubgraphs* reset(const CompressedSubgraphs& subgraphs, std::vector<VertexID>&& new_set) {
    keys_ = subgraphs.keys_;
    std::copy(subgraphs.sets_.begin(), subgraphs.sets_.end(), sets_.begin());
    sets_.back() = std::make_shared<std::vector<VertexID>>(std::move(new_set));
    return this;
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

  uint64_t getNumIsomorphicSubgraphs(const PruningIndexGroups& pruning_indices, uint64_t limit = ~0u) const;

  inline uint64_t getNumIsomorphicSubgraphs(uint64_t limit = ~0u) const {
    if (sets_.empty()) {
      return !keys_.empty();
    }
    std::vector<std::vector<VertexID>*> set_ptrs(sets_.size());
    for (uint32_t i = 0; i < sets_.size(); ++i) {
      set_ptrs[i] = sets_[i].get();
    }
    unordered_set<VertexID> existing_vertices;
    existing_vertices.reserve(getNumVertices() - 1);
    existing_vertices.insert(keys_.begin(), keys_.end());
    return getNumIsomorphicSubgraphs(existing_vertices, set_ptrs, limit);
  }

  /**
   * Count isomorphic subgraphs projected on the given sets
   */
  static uint64_t getNumIsomorphicSubgraphs(unordered_set<VertexID>& existing_vertices,
                                            std::vector<std::vector<VertexID>*>& set_ptrs, uint64_t limit = ~0u);

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

  std::ostream& logEnumerated(std::ostream& ss, const std::vector<std::pair<bool, uint32_t>>& log_indices,
                              uint64_t limit = ~0u) const;
  /** Get the value of the key vertex at key_idx. */
  VertexID getKeyVal(uint32_t key_idx) const {
    DCHECK_LT(key_idx, keys_.size()) << key_idx << "  " << keys_.size();
    return keys_[key_idx];
  }

  const std::vector<VertexID>& getKeys() const { return keys_; }

  std::vector<VertexID>& getKeys() { return keys_; }

  unordered_set<VertexID> getKeyMap() const { return unordered_set<VertexID>(keys_.begin(), keys_.end()); }

  // exceptions include keys and single-element sets with the same label
  inline void getExceptions(unordered_set<VertexID>& exception, const std::vector<uint32_t>& exception_key_indices,
                            const std::vector<uint32_t>& exception_set_indices) const {
    exception.reserve(exception_key_indices.size() + exception_set_indices.size());
    for (auto idx : exception_key_indices) {
      exception.insert(keys_[idx]);
    }

    for (auto idx : exception_set_indices) {
      DCHECK(sets_[idx] != nullptr) << idx << " " << keys_.size() << " + " << sets_.size();
      auto& set = *sets_[idx];
      if (set.size() == 1) {
        exception.insert(set.front());
      }
    }
  }

  // exceptions include keys and single-element sets with the same label
  inline unordered_set<VertexID> getExceptions(const std::vector<uint32_t>& exception_key_indices,
                                               const std::vector<uint32_t>& exception_set_indices) const {
    unordered_set<VertexID> exception;
    getExceptions(exception, exception_key_indices, exception_set_indices);
    return exception;
  }

  /** Get the matching set of the non-key vertex at key_idx. */
  const VertexSet& getSet(uint32_t key_idx) const { return sets_[key_idx]; }

  const std::vector<VertexSet>& getSets() const { return sets_; }

  void UpdateSets(uint32_t set_idx, VertexSet&& new_set) { sets_[set_idx] = std::move(new_set); }
  void UpdateSets(uint32_t set_idx, const VertexSet& new_set) { sets_[set_idx] = new_set; }

  /** Update the key vertex at key_idx to val. */
  void UpdateKey(uint32_t key_idx, VertexID val) { keys_[key_idx] = val; }
  /** Add a vertex val to the vertex set at set_idx. */
  void UpdateSet(uint32_t set_idx, VertexID val) { sets_[set_idx]->push_back(val); }

  bool empty() const { return keys_.empty(); }

  bool pruneExistingSets(VertexID v, unordered_set<uint32_t>& set_indices, uint32_t set_size_threshold,
                         bool recursive_prune = true);

  inline void swap(CompressedSubgraphs& other) {
    keys_.swap(other.keys_);
    sets_.swap(other.sets_);
  }
};

}  // namespace circinus

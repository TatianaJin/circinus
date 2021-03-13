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
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

/** A group of subgraphs compressed by some key vertices, which constitute a vertex cover for all subgraphs. */
class CompressedSubgraphs {
  std::vector<VertexID> keys_;
  std::vector<VertexSet> sets_;

 public:
  /**
   * @param key_size The number of key vertices.
   * @param n_vertices The number of vertices in each compressed subgraph.
   */
  [[deprecated]] CompressedSubgraphs(uint32_t key_size, uint32_t n_vertices)
      : keys_(key_size), sets_(n_vertices - key_size) {
    for (auto& set : sets_) {
      set = std::make_shared<std::vector<VertexID>>();
    }
  }

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
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, VertexID key)
      : keys_(subgraphs.getNumKeys() + 1), sets_(subgraphs.getNumSets()) {
    std::copy(subgraphs.keys_.begin(), subgraphs.keys_.end(), keys_.begin());
    keys_.back() = key;
    sets_ = subgraphs.sets_;
  }

  /**
   * @param subgraphs The compressed subgraphs that can extend to this CompressedSubraphs. They are one vertex smaller.
   * @param new_set The set of new vertices expanded in this CompressedSubgraphs, which is not in key.
   */
  CompressedSubgraphs(const CompressedSubgraphs& subgraphs, VertexSet&& new_set)
      : keys_(subgraphs.getNumKeys()), sets_(subgraphs.getNumSets() + 1) {
    keys_ = subgraphs.keys_;
    std::copy(subgraphs.sets_.begin(), subgraphs.sets_.end(), sets_.begin());
    sets_.back() = std::move(new_set);
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

  // TODO(tatiana): this function is quite expensive
  uint64_t getNumIsomorphicSubgraphs(uint64_t limit = ~0u) const {
    if (sets_.empty()) {
      return !keys_.empty();
    }
    if (sets_.size() == 1) {
      return sets_.front()->size();
    }
    // dfs sets_ chain
    uint64_t count = 0;
    std::vector<uint32_t> set_index(sets_.size(), 0);
    unordered_set<QueryVertexID> existing_vertices;
    existing_vertices.reserve(getNumVertices());
    existing_vertices.insert(keys_.begin(), keys_.end());
    // std::vector<QueryVertexID> existing_vertices;
    // existing_vertices.reserve(getNumVertices());
    // existing_vertices.insert(existing_vertices.end(), keys_.begin(), keys_.end());
    uint32_t last_depth = sets_.size() - 1;
    uint32_t current_depth = 0;
    while (true) {
      while (set_index[current_depth] < (*sets_[current_depth]).size()) {
        auto v = (*sets_[current_depth])[set_index[current_depth]];
        // bool existing = false;
        // for (auto e : existing_vertices) {
        //   if (e == v) {
        //     existing = true;
        //     break;
        //   }
        // }
        ++set_index[current_depth];
        // if (!existing) {                      // v is valid
        if (existing_vertices.count(v) == 0) {  // v is valid
          if (current_depth == last_depth) {    // reaching a leave in dfs
            if (++count == limit) return count;
          } else {
            // existing_vertices.push_back(v);
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
      // existing_vertices.pop_back();
      existing_vertices.erase((*sets_[current_depth])[set_index[current_depth] - 1]);
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

  /** Get the value of the key vertex at key_idx. */
  VertexID getKeyVal(uint32_t key_idx) const { return keys_[key_idx]; }

  unordered_set<VertexID> getKeyMap() const { return unordered_set<VertexID>(keys_.begin(), keys_.end()); }

  /** Get the matching set of the non-key vertex at key_idx. */
  const VertexSet& getSet(uint32_t key_idx) const { return sets_[key_idx]; }

  void UpdateSets(uint32_t set_idx, VertexSet&& new_set) { sets_[set_idx] = std::move(new_set); }

  /** Update the key vertex at key_idx to val. */
  void UpdateKey(uint32_t key_idx, VertexID val) { keys_[key_idx] = val; }
  /** Add a vertex val to the vertex set at set_idx. */
  void UpdateSet(uint32_t set_idx, VertexID val) { sets_[set_idx]->push_back(val); }
};

}  // namespace circinus

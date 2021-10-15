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

#include "graph/compressed_subgraphs.h"

#include <unordered_map>
#include <unordered_set>

#include "glog/logging.h"

namespace circinus {

bool CompressedSubgraphs::pruneExistingSets(VertexID v, unordered_set<uint32_t>& set_indices,
                                            uint32_t set_size_threshold, bool recursive_prune) {
  unordered_map<VertexID, uint32_t> new_v;  // vertex, set_index
  for (uint32_t i : set_indices) {
    auto& set = *sets_[i];
    if (set.size() <= set_size_threshold) {
      auto lb = circinus::lower_bound(set.begin(), set.end(), v);
      if (lb != set.end() && *lb == v) {  // conflict found
        if (set.size() == 1) {
          return true;
        }
        auto new_set = newSharedVSet();
        new_set->insert(new_set->end(), set.begin(), lb);
        new_set->insert(new_set->end(), lb + 1, set.end());
        // recursively prune when set size becomes 1
        if (new_set->size() == 1 && !new_v.insert({new_set->front(), i}).second) {
          return true;
        }
        sets_[i] = newVertexSet(new_set);
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

uint64_t CompressedSubgraphs::getNumIsomorphicSubgraphs(unordered_set<VertexID>& existing_vertices,
                                                        std::vector<const SingleRangeVertexSetView*>& set_ptrs,
                                                        uint64_t limit) {
  DCHECK(!set_ptrs.empty());
  std::sort(set_ptrs.begin(), set_ptrs.end(),
            [](const auto& set1, const auto& set2) { return set1->size() < set2->size(); });
  // dfs sets_ chain
  uint64_t count = 0;
  std::vector<uint32_t> set_index(set_ptrs.size(), 0);
  uint32_t last_depth = set_ptrs.size() - 1;
  uint32_t current_depth = 0;
  while (true) {
    while (set_index[current_depth] < (*set_ptrs[current_depth]).size()) {
      auto v = (*set_ptrs[current_depth])[set_index[current_depth]];
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
    existing_vertices.erase((*set_ptrs[current_depth])[set_index[current_depth] - 1]);
  }
  return count;
}
uint64_t CompressedSubgraphs::getNumIsomorphicSubgraphsWithConstraintsImpl(
    unordered_set<VertexID>& existing_vertices, std::vector<const SingleRangeVertexSetView*>& set_ptrs,
    const std::vector<std::vector<uint32_t>>& constraints_adj, uint64_t limit) {
  DCHECK(!set_ptrs.empty());
  // std::sort(set_ptrs.begin(), set_ptrs.end(),
  //           [](const auto& set1, const auto& set2) { return set1->size() < set2->size(); });
  // dfs sets_ chain
  uint64_t count = 0;
  std::vector<uint32_t> set_index(set_ptrs.size(), 0);
  uint32_t last_depth = set_ptrs.size() - 1;
  uint32_t current_depth = 0;
  while (true) {
    while (set_index[current_depth] < (*set_ptrs[current_depth]).size()) {
      auto v = (*set_ptrs[current_depth])[set_index[current_depth]];
      ++set_index[current_depth];
      if (existing_vertices.count(v) == 0) {  // v is valid
        if (current_depth == last_depth) {    // reaching a leave in dfs
          if (++count == limit) return count;
        } else {
          existing_vertices.insert(v);
          ++current_depth;
          uint64_t least_vid = 0;
          for (uint32_t nbr : constraints_adj[current_depth]) {
            least_vid = std::max(least_vid, (*set_ptrs[nbr])[set_index[nbr] - 1]);
          }
          set_index[current_depth] =
              std::lower_bound(set_ptrs[current_depth]->begin(), set_ptrs[current_depth]->end(), least_vid) -
              set_ptrs[current_depth]->begin();
        }
      }
    }
    if (current_depth == 0) {
      break;
    }
    --current_depth;
    existing_vertices.erase((*set_ptrs[current_depth])[set_index[current_depth] - 1]);
  }
  return count;
}

uint64_t CompressedSubgraphs::getNumIsomorphicSubgraphs(const PruningIndexGroups& pruning_indices,
                                                        uint64_t limit) const {
  if (pruning_indices.empty()) {
    return std::min(limit, getNumSubgraphs());
  }  // no same-label keys or sets
  if (sets_.empty()) {
    return !keys_.empty();
  }
  auto n_pruning_groups = pruning_indices.size();
  unordered_set<uint32_t> covered_indices;
  std::vector<std::vector<const SingleRangeVertexSetView*>> pruning_set_ptrs(n_pruning_groups);
  std::vector<uint32_t> pruning_set_group_size(n_pruning_groups, 1);
  for (uint32_t i = 0; i < n_pruning_groups; ++i) {
    auto pruning_set_indices = pruning_indices[i].second;
    DCHECK(!pruning_set_indices.empty());
    for (auto set_index : pruning_set_indices) {
      covered_indices.insert(set_index);
      auto& ptr = *getSet(set_index);
      pruning_set_ptrs[i].push_back(&ptr);
      pruning_set_group_size[i] *= ptr.size();
    }
  }

  uint64_t count = 1;
  // calculate count for sets not in the pruning groups
  if (getNumSets() > covered_indices.size()) {
    for (uint32_t i = 0; i < getNumSets(); ++i) {
      if (covered_indices.count(i) == 0) {
        // if (sets_[i]->empty()) return 0;
        count *= sets_[i]->size();
      }
    }
  }
  // calculate count for each pruning set group in ascending order of group size
  std::vector<uint32_t> pruning_set_group_index(n_pruning_groups);
  std::iota(pruning_set_group_index.begin(), pruning_set_group_index.end(), 0);
  std::sort(pruning_set_group_index.begin(), pruning_set_group_index.end(),
            [&pruning_set_group_size](uint32_t idx1, uint32_t idx2) {
              return pruning_set_group_size[idx1] < pruning_set_group_size[idx2];
            });

  unordered_set<VertexID> existing_vertices;
  for (auto idx : pruning_set_group_index) {
    existing_vertices.clear();
    for (auto key_i : pruning_indices[idx].first) {
      existing_vertices.insert(getKeyVal(key_i));
    }
    auto group_cnt = getNumIsomorphicSubgraphs(existing_vertices, pruning_set_ptrs[idx], (limit + count - 1) / count);
    if (group_cnt == 0) return 0;
    count *= group_cnt;
  }
  return std::min(count, limit);
}

uint64_t CompressedSubgraphs::getNumIsomorphicSubgraphsWithConstraints(
    const PruningIndexGroups& pruning_indices, const std::vector<std::vector<std::vector<uint32_t>>>& constraints_adjs,
    const std::vector<std::vector<uint32_t>>& enumerate_orders, uint64_t limit) const {
  // TODO(tatiana): reuse code in getNumIsomorphicSubgraphs >>>>>
  if (pruning_indices.empty()) {
    return std::min(limit, getNumSubgraphs());
  }  // no same-label keys or sets
  if (sets_.empty()) {
    return !keys_.empty();
  }
  auto n_pruning_groups = enumerate_orders.size();
  unordered_set<uint32_t> covered_indices;
  std::vector<std::vector<const SingleRangeVertexSetView*>> pruning_set_ptrs(n_pruning_groups);
  std::vector<uint32_t> pruning_set_group_size(n_pruning_groups, 1);
  for (uint32_t i = 0; i < n_pruning_groups; ++i) {
    DCHECK(!enumerate_orders[i].empty());
    for (auto set_index : enumerate_orders[i]) {
      covered_indices.insert(set_index);
      auto& ptr = *getSet(set_index);
      pruning_set_ptrs[i].push_back(&ptr);
      pruning_set_group_size[i] *= ptr.size();
    }
  }

  uint64_t count = 1;
  // calculate count for sets not in the pruning groups
  if (getNumSets() > covered_indices.size()) {
    for (uint32_t i = 0; i < getNumSets(); ++i) {
      if (covered_indices.count(i) == 0) {
        // if (sets_[i]->empty()) return 0;
        count *= sets_[i]->size();
      }
    }
  }
  // calculate count for each pruning set group in ascending order of group size
  std::vector<uint32_t> pruning_set_group_index(n_pruning_groups);
  std::iota(pruning_set_group_index.begin(), pruning_set_group_index.end(), 0);
  std::sort(pruning_set_group_index.begin(), pruning_set_group_index.end(),
            [&pruning_set_group_size](uint32_t idx1, uint32_t idx2) {
              return pruning_set_group_size[idx1] < pruning_set_group_size[idx2];
            });
  // TODO(tatiana): reuse code in getNumIsomorphicSubgraphs <<<<<

  unordered_set<VertexID> existing_vertices;
  for (auto idx : pruning_set_group_index) {
    existing_vertices.clear();
    for (auto key_i : pruning_indices[idx].first) {
      existing_vertices.insert(getKeyVal(key_i));
    }
    // FIXME(tatiana): use getNumIsomorphicSubgraphs if there is no constraint in this group
    auto group_cnt = getNumIsomorphicSubgraphsWithConstraintsImpl(existing_vertices, pruning_set_ptrs[idx],
                                                                  constraints_adjs[idx], (limit + count - 1) / count);
    if (group_cnt == 0) return 0;
    count *= group_cnt;
  }
  return std::min(count, limit);
}

std::ostream& CompressedSubgraphs::logEnumerated(std::ostream& ss,
                                                 const std::vector<std::pair<bool, uint32_t>>& log_indices,
                                                 uint64_t limit) const {
  if (sets_.empty()) {
    if (!keys_.empty()) {
      for (auto& pair : log_indices) {
        CHECK(pair.first);
        ss << keys_[pair.second] << ',';
      }
      ss << '\n';
    }
    return ss;
  }

  std::vector<VertexID> set_tuple(sets_.size(), -1);
  std::vector<std::vector<VertexID>*> set_ptrs(sets_.size());
  for (uint32_t i = 0; i < sets_.size(); ++i) {
    set_ptrs[i] = sets_[i].get();
  }
  unordered_set<VertexID> existing_vertices;
  existing_vertices.reserve(getNumVertices() - 1);
  existing_vertices.insert(keys_.begin(), keys_.end());
  std::vector<uint32_t> set_indices(sets_.size());
  std::iota(set_indices.begin(), set_indices.end(), 0);
  std::sort(set_indices.begin(), set_indices.end(),
            [&set_ptrs](uint32_t set1, uint32_t set2) { return set_ptrs[set1]->size() < set_ptrs[set2]->size(); });
  std::sort(set_ptrs.begin(), set_ptrs.end(),
            [](const auto& set1, const auto& set2) { return set1->size() < set2->size(); });
  // dfs sets_ chain
  uint64_t count = 0;
  std::vector<uint32_t> set_index(set_ptrs.size(), 0);
  uint32_t last_depth = set_ptrs.size() - 1;
  uint32_t current_depth = 0;
  while (true) {
    while (set_index[current_depth] < (*set_ptrs[current_depth]).size()) {
      auto v = (*set_ptrs[current_depth])[set_index[current_depth]];
      ++set_index[current_depth];
      if (existing_vertices.count(v) == 0) {  // v is valid
        set_tuple[set_indices[current_depth]] = v;
        if (current_depth == last_depth) {  // reaching a leave in dfs
          for (auto& pair : log_indices) {
            if (pair.first) {
              ss << keys_[pair.second] << ',';
            } else {
              ss << set_tuple[pair.second] << ',';
            }
          }
          ss << '\n';
          if (++count == limit) return ss;
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
    existing_vertices.erase((*set_ptrs[current_depth])[set_index[current_depth] - 1]);
  }
  return ss;
}

}  // namespace circinus

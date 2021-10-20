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

#include "algorithms/intersect.h"

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

uint64_t CompressedSubgraphs::getNumIsomorphicSubgraphsImpl(unordered_set<VertexID>& existing_vertices,
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

// TODO(tatiana): compute an inclusion-exclusion schedule with intersection and coefficient and pass the schedule to
// this function
uint64_t CompressedSubgraphs::getNumIsomorphicSubgraphs(const PruningIndexGroups& pruning_indices, uint64_t limit,
                                                        const VertexRelationship* qv_equivalent_classes) const {
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
  VertexID max_existing_vertices = 0;
  unordered_map<std::string, std::vector<VertexID>> cache;
  for (auto idx : pruning_set_group_index) {
    existing_vertices.clear();
    max_existing_vertices = 0;
    for (auto key_i : pruning_indices[idx].first) {
      existing_vertices.insert(getKeyVal(key_i));
      max_existing_vertices = std::max(max_existing_vertices, getKeyVal(key_i));
    }
    uint64_t group_cnt = 0;
    if (pruning_set_ptrs[idx].size() == 1) {  // existing_vertices must not be empty
      group_cnt = countExcept(*pruning_set_ptrs[idx].front(), existing_vertices, max_existing_vertices);
    } else if (pruning_set_ptrs[idx].size() < 4) {
      group_cnt = getEnumerationCount(pruning_indices[idx].second, qv_equivalent_classes, cache, existing_vertices,
                                      max_existing_vertices);
    } else {
      group_cnt = getNumIsomorphicSubgraphsImpl(existing_vertices, pruning_set_ptrs[idx], (limit + count - 1) / count);
    }
    if (group_cnt == 0) return 0;
    count *= group_cnt;
    cache.clear();
  }
  return std::min(count, limit);
}

uint64_t CompressedSubgraphs::countExcept(const SingleRangeVertexSetView& set, const unordered_set<VertexID>& except,
                                          VertexID max_except) {
  uint64_t count = 0;
  for (auto it = set.begin(); it != set.end(); ++it) {
    count += (except.count(*it) == 0);
    if (*it > max_except) {  // all the following are not excepted
      return count + std::distance(++it, set.end());
    }
  }
  return count;
}

uint64_t CompressedSubgraphs::getEnumerationCount(const std::vector<uint32_t>& set_indices,
                                                  const VertexRelationship* qv_equivalent_classes,
                                                  unordered_map<std::string, std::vector<VertexID>>& cache,
                                                  const unordered_set<VertexID>& except, VertexID max_except) const {
  uint32_t n_sets = set_indices.size();
  DCHECK_GT(n_sets, 1) << set_indices[0];
  DCHECK_LT(n_sets, 4) << "counting more than 3 sets is not supported yet";  // TODO(tatiana)
  uint64_t cnt = 0, product = 1;
  if (except.empty()) {
    for (auto index : set_indices) {
      product *= getSet(index)->size();
    }
    cnt = product;
    uint64_t size = 0;
    for (uint32_t i = 0; i < n_sets; ++i) {
      for (uint32_t j = i + 1; j < n_sets; ++j) {
        bool to_cache = (j + 1) < n_sets;
        size = getIntersectionSize(set_indices[i], set_indices[j], qv_equivalent_classes, to_cache ? &cache : nullptr,
                                   except);
        size *= product / getSet(set_indices[i])->size() / getSet(set_indices[j])->size();
        cnt -= size;
      }
    }
  } else {
    // actually no need to use except if the union of neighbors of sets covers all same label keys
    std::vector<uint32_t> set_sizes;
    set_sizes.reserve(n_sets);
    for (auto index : set_indices) {
      set_sizes.push_back(countExcept(*getSet(index), except, max_except));
      product *= set_sizes.back();
    }
    cnt = product;
    for (uint32_t i = 0; i < n_sets; ++i) {
      for (uint32_t j = i + 1; j < n_sets; ++j) {
        bool to_cache = (j + 1) < n_sets;
        auto size = getIntersectionSize(set_indices[i], set_indices[j], qv_equivalent_classes,
                                        to_cache ? &cache : nullptr, except);
        size *= product / set_sizes[set_indices[i]] / set_sizes[set_indices[j]];
        cnt -= size;
      }
    }
    /*
    CHECK_EQ(truth, cnt) << "max_existing_vertices " << max_except << " sets " << set_sizes[set_indices[0]] << '/'
                         << getSet(set_indices[0])->size() << ' ' << set_sizes[set_indices[1]] << '/'
                         << getSet(set_indices[1])->size() << "product " << product << " intersection "
                         << getIntersectionSize(set_indices[0], set_indices[1], qv_equivalent_classes, nullptr, except)
                         << " coeff " << product / set_sizes[set_indices[0]] / set_sizes[set_indices[1]];
                         */
  }
  if (n_sets == 3) {
    cnt += 2 * intersectionCount(cache.begin()->second, *getSet(set_indices[2]));  // no need to except again
  }
  return cnt;
}

uint64_t CompressedSubgraphs::getIntersectionSize(uint32_t set_idx1, uint32_t set_idx2,
                                                  const VertexRelationship* qv_equivalent_classes,
                                                  unordered_map<std::string, std::vector<VertexID>>* cache,
                                                  const unordered_set<VertexID>& except) const {
  if (qv_equivalent_classes->isEquivalentByIndices(set_idx1, set_idx2)) {
    auto n = getSet(set_idx1)->size();
    DCHECK_EQ(n, getSet(set_idx2)->size()) << "set indices" << set_idx1 << ' ' << set_idx2;
    return n;
  }
  // cross product - intersection
  auto& set1 = *getSet(set_idx1);
  auto& set2 = *getSet(set_idx2);
  if (cache == nullptr) {
    return intersectionCount(set1, set2, except);
  }
  std::vector<VertexID> intersection;
  intersect(set1, set2, &intersection, except);
  auto size = intersection.size();
  // cache the result
  if (set_idx1 > set_idx2) {
    std::swap(set_idx1, set_idx2);
  }
  uint64_t key = (set_idx1 << sizeof(uint32_t)) | set_idx2;
  cache->insert({std::string((char*)&key, sizeof(key)), std::move(intersection)});
  return size;
}

uint64_t CompressedSubgraphs::getNumIsomorphicSubgraphsWithConstraints(
    const PruningIndexGroups& pruning_indices, const std::vector<std::vector<std::vector<uint32_t>>>& constraints_adjs,
    const std::vector<std::vector<uint32_t>>& enumerate_orders, uint64_t limit,
    const VertexRelationship* qv_equivalent_classes) const {
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
  VertexID max_existing_vertices = 0;
  unordered_map<std::string, std::vector<VertexID>> cache;
  for (auto idx : pruning_set_group_index) {
    existing_vertices.clear();
    max_existing_vertices = 0;
    for (auto key_i : pruning_indices[idx].first) {
      existing_vertices.insert(getKeyVal(key_i));
      max_existing_vertices = std::max(max_existing_vertices, getKeyVal(key_i));
    }
    uint64_t group_cnt;
    if (pruning_set_ptrs[idx].size() == 1) {  // existing_vertices must not be empty
      group_cnt = countExcept(*pruning_set_ptrs[idx].front(), existing_vertices, max_existing_vertices);
    } else if (pruning_set_ptrs[idx].size() == 2 &&
               qv_equivalent_classes->isEquivalentByIndices(enumerate_orders[idx].front(),
                                                            enumerate_orders[idx].back())) {
      auto size = pruning_set_ptrs[idx].front()->size();
      DCHECK_EQ(size, pruning_set_ptrs[idx].back()->size());  // two sets shall be the same
      group_cnt = size * (size - 1);
      if (!constraints_adjs[idx].back().empty()) {  // partial order between them
        group_cnt /= 2;
      }
    } else if (pruning_set_ptrs[idx].size() < 4) {  // sets without partial order constraint
      bool use_inex = true;
      for (auto& conds : constraints_adjs[idx]) {
        if (!conds.empty()) {
          use_inex = false;
          break;
        }
      }
      if (use_inex) {
        group_cnt = getEnumerationCount(enumerate_orders[idx], qv_equivalent_classes, cache, existing_vertices,
                                        max_existing_vertices);
      } else {
        group_cnt = getNumIsomorphicSubgraphsWithConstraintsImpl(existing_vertices, pruning_set_ptrs[idx],
                                                                 constraints_adjs[idx], (limit + count - 1) / count);
      }
    } else {
      // FIXME(tatiana): use getNumIsomorphicSubgraphs if there is no constraint in this group
      group_cnt = getNumIsomorphicSubgraphsWithConstraintsImpl(existing_vertices, pruning_set_ptrs[idx],
                                                               constraints_adjs[idx], (limit + count - 1) / count);
    }
    if (group_cnt == 0) return 0;
    count *= group_cnt;
    cache.clear();
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

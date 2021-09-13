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
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/search.h"
#include "graph/types.h"
#include "graph/vertex_set_view.h"
#include "utils/hashmap.h"

namespace circinus {

template <typename Set1, typename Set2, typename V = VertexID>
inline void intersectIter(const std::enable_if_t<!std::is_same_v<Set1, VertexSetView>, Set1>& set1, const Set2& set2,
                          std::vector<V>* intersection, const unordered_set<V>& except) {
  // LOG(FATAL) << "should not enter this function " << getTypename(set1) << ", " << getTypename(set2);
  intersection->reserve(set1.size());
  if
    constexpr(std::is_same_v<Set2, VertexSetView>) {
      if (set2.getRanges().size() == 1) {
        auto lower_bound = set2.getRanges().front().first;
        auto end = lower_bound + set2.getRanges().front().second;
        for (V vid : set1) {
          lower_bound = lowerBound(lower_bound, end, vid);
          // if lower_bound is out of range, all elements in the rest of set2 are smaller than the rest of set1
          if (lower_bound == end) {
            break;
          }
          if (*lower_bound == vid && except.count(vid) == 0) intersection->emplace_back(vid);
        }
        return;
      }
    }

  auto lower_bound = set2.begin();
  for (V vid : set1) {
    lower_bound = lowerBound(lower_bound, set2.end(), vid);
    // if lower_bound is out of range, all elements in the rest of set2 are smaller than the rest of set1
    if (lower_bound == set2.end()) {
      break;
    }
    if (*lower_bound == vid && except.count(vid) == 0) intersection->emplace_back(vid);
  }
}

template <typename Set1, typename Set2, typename V = VertexID>
inline void intersectIter(const std::enable_if_t<std::is_same_v<Set1, VertexSetView>, Set1>& set1, const Set2& set2,
                          std::vector<V>* intersection, const unordered_set<V>& except) {
  intersection->reserve(set1.size());
  if
    constexpr(std::is_same_v<Set2, VertexSetView>) {
      if (set2.getRanges().size() == 1) {
        auto lower_bound = set2.getRanges().front().first;
        auto end = lower_bound + set2.getRanges().front().second;
        for (auto& range : set1.getRanges()) {
          for (uint32_t i = 0; i < range.second; ++i) {
            auto vid = range.first[i];
            lower_bound = lowerBound(lower_bound, end, vid);
            // if lower_bound is out of range, all elements in the rest of set2 are smaller than the rest of set1
            if (lower_bound == end) {
              break;
            }
            if (*lower_bound == vid && except.count(vid) == 0) intersection->emplace_back(vid);
          }
        }
        return;
      }
    }

  auto lower_bound = set2.begin();
  for (auto& range : set1.getRanges()) {
    for (uint32_t i = 0; i < range.second; ++i) {
      auto vid = range.first[i];
      lower_bound = lowerBound(lower_bound, set2.end(), vid);
      // if lower_bound is out of range, all elements in the rest of set2 are smaller than the rest of set1
      if (lower_bound == set2.end()) {
        break;
      }
      if (*lower_bound == vid && except.count(vid) == 0) intersection->emplace_back(vid);
    }
  }
}

/** The given set1 and set2 must be sorted in ascending order.
 *
 * Set1 and Set2 should be forward iterable and provide the size() function.
 */
template <typename Set1, typename Set2, typename V = VertexID>
inline void intersect(const Set1& set1, const Set2& set2, std::vector<V>* intersection,
                      const unordered_set<V>& except = {}) {
  if (set1.size() > set2.size()) {
    return intersectIter<Set2, Set1>(set2, set1, intersection, except);
  }
  return intersectIter<Set1, Set2>(set1, set2, intersection, except);
}

/**
 * Set should be forward iterable and provide the size() function
 */
template <typename Set, typename V = VertexID>
inline void intersect(const unordered_set<V>& set1, const Set& set2, std::vector<V>* intersection,
                      const unordered_set<V>& except = {}) {
  intersection->reserve(std::min((size_t)set2.size(), set1.size()));
  for (auto vid : set2) {
    if (set1.count(vid) && except.count(vid) == 0) intersection->emplace_back(vid);
  }
}

/**
 * The given set1 and set2 must be sorted in ascending order.
 * Set should be forward iterable and provide the size() function
 *
 * @param intersection The non-const pointer to set1.
 */
template <typename Set, typename V = VertexID>
void intersectInplace(const std::vector<V>& set1, const Set& set2, std::vector<V>* intersection) {
  uint32_t size = 0;
  auto lb = set2.begin();
  auto end = set2.end();
  for (auto vid : set1) {
    lb = lowerBound(lb, end, vid);
    if (lb == end) {
      break;
    }
    if (*lb == vid) {
      (*intersection)[size++] = vid;
    }
  }
  intersection->resize(size);
}

inline void intersectInplace(std::vector<VertexID>* set1, const unordered_set<VertexID>& set2) {
  uint32_t size = 0;
  for (auto vid : *set1) {
    if (set2.count(vid)) {
      (*set1)[size++] = vid;
    }
  }
  set1->resize(size);
}

}  // namespace circinus
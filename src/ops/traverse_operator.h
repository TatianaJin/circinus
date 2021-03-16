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
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/types.h"
#include "ops/operator.h"
#include "utils/hashmap.h"

namespace circinus {

/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::pair<const VertexID*, uint32_t>& set1,
                      const std::pair<const VertexID*, uint32_t>& set2, std::vector<VertexID>* intersection,
                      const unordered_set<VertexID>& except = {}) {
  if (set1.second <= set2.second) {
    intersection->reserve(set1.second);
    uint32_t set2_index = 0;
    for (uint32_t i = 0; i < set1.second; ++i) {
      auto vid = set1.first[i];
      while (set2_index < set2.second && set2.first[set2_index] < vid) {
        ++set2_index;
      }
      if (set2_index == set2.second) {  // all elements in the rest of set2 are smaller than the rest of set1
        break;
      }
      if (set2.first[set2_index] == vid && except.count(vid) == 0) intersection->emplace_back(vid);
    }
  } else {
    intersect(set2, set1, intersection, except);
  }
}

inline void intersect(const unordered_set<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection, const unordered_set<VertexID>& except = {}) {
  intersection->reserve(std::min(set2.second, (uint32_t)set1.size()));
  for (uint32_t i = 0; i < set2.second; ++i) {
    auto vid = set2.first[i];
    if (set1.count(vid) && except.count(vid) == 0) intersection->emplace_back(vid);
  }
}

/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection, const unordered_set<VertexID>& except = {}) {
  intersect(std::make_pair(set1.data(), (uint32_t)set1.size()), set2, intersection, except);
}

/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::vector<VertexID>& set1, const std::vector<VertexID>& set2,
                      std::vector<VertexID>* intersection, const unordered_set<VertexID>& except = {}) {
  intersect(set1, std::make_pair(set2.data(), (uint32_t)set2.size()), intersection, except);
}

void intersectInplace(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection);

class TraverseOperator : public Operator {
 protected:
  const std::vector<VertexID>* candidates_ = nullptr;

  /* transient variables for recording the current inputs */
  uint32_t input_index_ = 0;
  const std::vector<CompressedSubgraphs>* current_inputs_ = nullptr;
  const Graph* current_data_graph_ = nullptr;

 public:
  virtual ~TraverseOperator() {}

  inline virtual void setCandidateSets(const std::vector<VertexID>* candidates) { candidates_ = candidates; }
  inline const std::vector<VertexID>* getCandidateSets() const { return candidates_; }

  virtual void input(const std::vector<CompressedSubgraphs>& inputs, const Graph* data_graph) {
    current_inputs_ = &inputs;
    input_index_ = 0;
    current_data_graph_ = data_graph;
  }

  virtual uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) = 0;
};

}  // namespace circinus

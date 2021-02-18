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

#include "ops/traverse_operator.h"

#include <utility>
#include <vector>

#include "graph/types.h"

namespace circinus {

void intersect(const std::pair<const VertexID*, uint32_t>& set1, const std::pair<const VertexID*, uint32_t>& set2,
               std::vector<VertexID>* intersection) {
  if (set1.second < set2.second) {
    auto lower_bound = set2.first;
    for (uint32_t i = 0; i < set1.second; ++i) {
      auto vid = set1.first[i];
      lower_bound = std::lower_bound(lower_bound, set2.first + set2.second, vid);
      uint32_t index = lower_bound - set2.first;
      // if lower_bound is out of range, all elements in the rest of set2 are smaller than the rest of set1
      if (index >= set2.second) {
        break;
      }
      if (*lower_bound == vid) intersection->emplace_back(vid);
    }
  } else {
    intersect(set2, set1, intersection);
  }
}

}  // namespace circinus

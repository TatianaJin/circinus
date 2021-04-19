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

#include "ops/traverse_operator.h"

#include <utility>
#include <vector>

#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

void intersectInplace(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection) {
  uint32_t size = 0;
  uint32_t lb = 0;
  for (uint32_t i = 0; i < set1.size(); ++i) {
    VertexID vid = set1[i];
    while (lb < set2.second && set2.first[lb] < vid) {
      ++lb;
    }
    if (lb == set2.second) {
      break;
    }

    if (set2.first[lb] == vid) {
      (*intersection)[size++] = vid;
    }
  }
  intersection->resize(size);
}

void intersectInplace_bs(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                         std::vector<VertexID>* intersection) {
  uint32_t size = 0;
  uint32_t lb = 0;
  for (uint32_t i = 0; i < set1.size(); ++i) {
    VertexID vid = set1[i];
    lb = std::lower_bound(set2.first + lb, set2.first + set2.second, vid) - set2.first;

    if (lb == set2.second) {
      break;
    }

    if (set2.first[lb] == vid) {
      (*intersection)[size++] = vid;
    }
  }
  intersection->erase(intersection->begin() + size, intersection->end());
}

}  // namespace circinus

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

#include <vector>

#include "graph/graph.h"
#include "graph/types.h"

namespace circinus {

class LocalFilter {
 public:
  virtual ~LocalFilter() {}

  /** @returns The number of records that passed the filter and are added to output */
  uint32_t filter(const Graph& data_graph, const std::vector<VertexID>& candidates,
                  std::vector<VertexID>* output) const {
    uint32_t count = 0;
    for (auto candidate : candidates) {
      if (prune(data_graph, candidate)) {
        continue;
      }
      ++count;
      output->push_back(candidate);
    }
    return count;
  }

  /** @returns True if vertex v is not a valid mapping */
  virtual bool prune(const Graph& data_graph, VertexID v) const = 0;
};

}  // namespace circinus

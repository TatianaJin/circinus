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

namespace circinus {

class ExpandVertexTraverseContext : public TraverseContext {
  // for profiling
  std::vector<unordered_set<std::string>> parent_tuple_sets_;

 public:
  explicit ExpandVertexTraverseContext(uint32_t input_index, uint32_t input_end_index,
                                       const std::vector<CompressedSubgraphs>* inputs, const void* data_graph,
                                       uint32_t parent_size)
      : TraverseContext(input_index, input_end_index, inputs, data_graph), parent_tuple_sets_(parent_size) {}

  inline void updateDistinctSICount(uint32_t depth, std::vector<VertexID>& parent_tuple, uint32_t pidx) {
    distinct_intersection_count +=
        parent_tuple_sets_[depth].emplace((const char*)parent_tuple.data(), (pidx + 1) * sizeof(VertexID)).second;
  }
};

}  // namespace circinus
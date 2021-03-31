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

#include "ops/expand_set_to_key_vertex_operator.h"

namespace circinus {

void ExpandSetToKeyVertexOperator::updateDistinctSICount() {
  auto& input = (*current_inputs_)[input_index_];
  std::vector<std::vector<VertexID>*> parent_set_ptrs;
  parent_set_ptrs.reserve(parents_.size());
  for (auto parent : parents_) {
    parent_set_ptrs.push_back(input.getSet(query_vertex_indices_[parent]).get());
  }
  uint32_t depth = 0, last_depth = parents_.size() - 1;
  std::vector<uint32_t> set_index(parents_.size(), 0);
  std::vector<VertexID> parent_tuple(parents_.size());
  unordered_set<VertexID> prefix_set;
  while (true) {
    while (set_index[depth] < parent_set_ptrs[depth]->size()) {
      auto parent_vid = (*parent_set_ptrs[depth])[set_index[depth]];
      if (prefix_set.count(parent_vid)) {  // the parent tuples with current prefix will be pruned
        ++set_index[depth];
        continue;
      }
      parent_tuple[depth] = parent_vid;
      distinct_intersection_count_ +=
          parent_tuple_sets_[depth].emplace((char*)parent_tuple.data(), (depth + 1) * sizeof(VertexID)).second;
      if (depth == last_depth) {
        ++set_index[depth];
      } else {
        prefix_set.insert(parent_vid);
        ++depth;
        set_index[depth] = 0;
      }
    }
    if (depth == 0) {
      break;
    }
    --depth;
    prefix_set.erase((*parent_set_ptrs[depth])[set_index[depth]]);
    ++set_index[depth];
  }
}

}  // namespace circinus

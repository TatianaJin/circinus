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

#include <string>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/query_graph.h"
#include "ops/traverse_operator.h"

namespace circinus {

class ExpandVertexOperator : public TraverseOperator {
 public:
  ExpandVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                       const std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices)
      : parents_(parents), target_vertex_(target_vertex), query_vertex_indices_(query_vertex_indices) {}

 protected:
  inline void toStringInner(std::stringstream& ss) const {
    for (auto parent : parents_) {
      DCHECK_EQ(query_vertex_indices_.count(parent), 1);
      ss << ' ' << parent;
    }
    DCHECK_EQ(query_vertex_indices_.count(target_vertex_), 1);
    ss << " -> " << target_vertex_;
    if (candidates_ != nullptr) ss << " (" << candidates_->size() << ")";
  }

  std::vector<QueryVertexID> parents_;
  std::unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;
  QueryVertexID target_vertex_;
};

}  // namespace circinus

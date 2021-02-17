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

#include "graph/query_graph.h"
#include "ops/traverse_operator.h"

namespace circinus {

class ExpandVertexOperator : public TraverseOperator {
 public:
  ExpandVertexOperator(const QueryGraph* g, std::vector<QueryVertexID>& parent_vertex, QueryVertexID target_vertex,
                       const std::vector<int>& cover_table) {
    // FIXME
  }
};

}  // namespace circinus

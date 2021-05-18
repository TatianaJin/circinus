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

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters/filter.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

namespace circinus {

class GQLFilter : public Filter {
 private:
  bool verify(const Graph& data_graph, VertexID data_vertex);

 public:
  GQLFilter(const QueryGraph* query_graph, const Graph* data_graph, QueryVertexID query_vertex);

  void filter(const Graph& data_graph, std::vector<VertexID>& candidates);
};

}  // namespace circinus

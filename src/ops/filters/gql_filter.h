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

class GQLFilter : public NeighborhoodFilter {
 private:
  bool verify(const GraphBase* data_graph, const VertexID data_vertex,
              std::vector<std::vector<VertexID>>* candidates) const;

 public:
  GQLFilter(ExecutionConfig& conf, const QueryGraph* query_graph, QueryVertexID query_vertex)
      : NeighborhoodFilter(conf, query_graph, query_vertex) {}

  void filter(const GraphBase* data_graph, std::vector<std::vector<VertexID>>* candidates,
              FilterContext* ctx) const override;
};

}  // namespace circinus

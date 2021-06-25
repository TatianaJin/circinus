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
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "ops/order/order_base.h"
#include "algorithms/k_core.h"

namespace circinus {

class CFLOrder : public OrderBase {
 public:
  TwoCoreSolver two_core_solver_;
  QueryVertexID getStartVertex(const GraphMetadata& metadata, const QueryGraph* query_graph,
                               const std::vector<VertexID>& candidate_size) override;

  std::vector<QueryVertexID> getTopThree(const GraphMetadata& metadata, const QueryGraph* q);

  QueryVertexID getStartVertex(const std::vector<QueryVertexID>& query_vertices,
                               const std::vector<VertexID>& cardinality, const QueryGraph& q,
                               const GraphMetadata& metadata);
};

}  // namespace circinus

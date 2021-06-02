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

#include "ops/logical/filter/gql_filter.h"

#include <memory>

#include "ops/filters/filter.h"
#include "ops/filters/gql_filter.h"
#include "utils/utils.h"

namespace circinus {

LogicalGQLFilter::LogicalGQLFilter(const QueryGraph* query_graph) : LogicalNeighborhoodFilter(query_graph) {}

std::vector<std::unique_ptr<NeighborhoodFilter>> LogicalGQLFilter::toPhysicalOperators(const GraphMetadata& metadata,
                                                                                       ExecutionConfig& exec) {
  std::vector<std::unique_ptr<NeighborhoodFilter>> ret;
  for (uint32_t i = 0; i < 2; ++i) {
    for (QueryVertexID query_vertex = 0; query_vertex < query_graph_->getNumVertices(); ++query_vertex) {
      ret.emplace_back(std::make_unique<GQLFilter>(exec, query_graph_, query_vertex));
    }
  }
  return ret;
}

}  // namespace circinus

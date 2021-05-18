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

#include <algorithm>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "utils/hashmap.h"

namespace circinus {

class Filter;

class LogicalFilter {
 protected:
  const QueryGraph* query_graph_;
  const Graph* data_graph_;

 public:
  LogicalFilter(const QueryGraph* query_graph, const Graph* data_graph)
      : query_graph_(query_graph), data_graph_(data_graph) {}

  virtual ~LogicalFilter() {}

  virtual std::vector<std::unique_ptr<Filter>> toPhysicalOperators(const GrapMetadata& metadata,
                                                                   ExecutionConfig& exec) = 0;
};

}  // namespace circinus

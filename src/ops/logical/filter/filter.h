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

#include <memory>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph_metadata.h"
#include "graph/types.h"

namespace circinus {

class LocalFilter;         // forward declaration
class NeighborhoodFilter;  // forward declaration
class Extender;            // forward declaration

class LogicalLocalFilter {
 public:
  virtual ~LogicalLocalFilter() {}
  virtual std::vector<std::unique_ptr<LocalFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                        ExecutionConfig& exec) = 0;
};

class LogicalNeighborhoodFilter {
 protected:
  const QueryGraph* query_graph_;

 public:
  explicit LogicalNeighborhoodFilter(const QueryGraph* query_graph) : query_graph_(query_graph) {}

  virtual ~LogicalNeighborhoodFilter() {}
  virtual std::vector<std::unique_ptr<NeighborhoodFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                               ExecutionConfig& exec) = 0;
  // virtual std::vector<std::unique_ptr<Extender>> toPhysicalExtenders(const GraphMetadata& metadata,
  //                                                                    ExecutionConfig& exec) = 0;
};
}  // namespace circinus

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
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

#include "glog/logging.h"

namespace circinus {

// forward declaration
class InputOperator;

class LogicalCompressedInputOperator {
 protected:
  const QueryVertexID input_query_vertex_;
  const bool inputs_are_keys_;

 public:
  LogicalCompressedInputOperator(const QueryVertexID input_query_vertex, const bool inputs_are_keys)
      : input_query_vertex_(input_query_vertex), inputs_are_keys_(inputs_are_keys) {}

  virtual ~LogicalCompressedInputOperator() {}

  virtual std::unique_ptr<InputOperator> toPhysicalOperators();
};

class PartitionedLogicalCompressedInputOperator : public LogicalCompressedInputOperator {
 protected:
  std::vector<std::pair<QueryVertexID, QueryVertexID>> qv_pivots_;

 public:
  PartitionedLogicalCompressedInputOperator(const QueryGraph* query_graph, bool inputs_are_keys,
                                            const std::vector<QueryVertexID>& matching_order,
                                            const std::vector<QueryVertexID>& partitioning_qvs);

  virtual ~PartitionedLogicalCompressedInputOperator() {}

  std::unique_ptr<InputOperator> toPhysicalOperators() override;
};

}  // namespace circinus

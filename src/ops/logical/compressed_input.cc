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

#include "ops/logical/compressed_input.h"

#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/input_operator.h"
#include "utils/hashmap.h"

#include "glog/logging.h"

namespace circinus {

LogicalCompressedInputOperator::LogicalCompressedInputOperator(bool inputs_are_keys,
                                                               const std::vector<QueryVertexID>& matching_order,
                                                               const std::vector<QueryVertexID>& partitioning_qvs)
    : input_query_vertex_(matching_order.front()), inputs_are_keys_(inputs_are_keys) {
  // find the last partitioning query vertex in order
  unordered_set<QueryVertexID> partitioning_qv_set(partitioning_qvs.begin(), partitioning_qvs.end());
  uint32_t pivot_order = matching_order.size() - 1;
  for (; pivot_order > 0; --pivot_order) {
    if (partitioning_qv_set.count(matching_order[pivot_order]) == 1) {
      break;
    }
  }
  DCHECK_EQ(partitioning_qv_set.count(matching_order[pivot_order]), 1);

  // if there is only one partitioning query vertex, and it is the input query vertex, no pruning is needed
  // otherwise, prune candidates for the first query vertex by reverse matching order
  if (pivot_order == 0) {
    return;
  }
  qv_pivots_.reserve(pivot_order);
  for (; pivot_order > 0; --pivot_order) {
    qv_pivots_.push_back({matching_order[pivot_order - 1], {matching_order[pivot_order]}});
  }
}

std::vector<std::unique_ptr<InputOperator>> LogicalCompressedInputOperator::toPhysicalOperators() {
  return {};  // TODO(byli)
}

}  // namespace circinus


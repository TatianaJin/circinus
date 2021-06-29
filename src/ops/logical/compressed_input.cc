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

#include "ops/logical/compressed_input.h"

#include <utility>
#include <vector>

#include "graph/tree_node.h"
#include "graph/types.h"
#include "ops/input_operator.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

#include "glog/logging.h"

namespace circinus {

LogicalCompressedInputOperator::LogicalCompressedInputOperator(const QueryGraph* query_graph, bool inputs_are_keys,
                                                               const std::vector<QueryVertexID>& matching_order,
                                                               const std::vector<QueryVertexID>& partitioning_qvs)
    : input_query_vertex_(matching_order.front()), inputs_are_keys_(inputs_are_keys) {
  // find the last partitioning query vertex in order
  unordered_set<QueryVertexID> partitioning_qv_set(partitioning_qvs.begin(), partitioning_qvs.end());

  if (partitioning_qv_set.count(input_query_vertex_) == 1) {
    return;
  }

  std::queue<QueryVertexID> bfs_queue;
  std::vector<bool> visited(matching_order.size(), false);
  std::vector<QueryVertexID> father;
  std::vector<QueryVertexID> distance(matching_order.size(), 0);

  bfs_queue.push(input_query_vertex_);
  visited[input_query_vertex_] = true;
  QueryVertexID pv = 0;
  // bfs to find the path from a partitioning vertex with minimal level
  while (!bfs_queue.empty()) {
    QueryVertexID u = bfs_queue.front();
    bfs_queue.pop();
    auto nbrs = query_graph->getOutNeighbors(u);
    bool find_pv = false;
    for (uint32_t i = 0; i < nbrs.second; ++i) {
      QueryVertexID v = nbrs.first[i];
      if (visited[v] == true) {
        continue;
      }
      father[v] = u;
      if (partitioning_qv_set.count(v) == 1) {
        find_pv = true;
        pv = v;
        break;
      }
      visited[v] = true;
      bfs_queue.push(v);
    }
    if (find_pv) {
      break;
    }
  }
  qv_pivots_.reserve(matching_order.size());
  while (true) {
    qv_pivots_.push_back({father[pv], pv});
    pv = father[pv];
    if (pv == input_query_vertex_) {
      break;
    }
  }
}

std::unique_ptr<InputOperator> LogicalCompressedInputOperator::toPhysicalOperators() {
  return std::make_unique<PartitionedInputOperator>(input_query_vertex_, inputs_are_keys_, &qv_pivots_);
}

}  // namespace circinus

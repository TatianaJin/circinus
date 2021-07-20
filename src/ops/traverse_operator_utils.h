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

#include <utility>

#include "graph/bipartite_graph.h"
#include "graph/graph.h"
#include "graph/graph_view.h"
#include "graph/partitioned_graph.h"
#include "graph/types.h"
#include "ops/traverse_operator.h"

namespace circinus {

template <template <typename> typename OpType, typename... Args>
inline TraverseOperator* newTraverseOp(GraphType g_type, Args... args) {
  switch (g_type) {
  case GraphType::Normal:
    return new OpType<Graph>(std::forward<Args>(args)...);
  case GraphType::Partitioned:
    return new OpType<ReorderedPartitionedGraph>(std::forward<Args>(args)...);
  case GraphType::GraphView:
    return new OpType<GraphView<GraphPartitionBase>>(std::forward<Args>(args)...);
  case GraphType::BipartiteGraphView:
    return new OpType<GraphView<BipartiteGraph>>(std::forward<Args>(args)...);
  }
  LOG(FATAL) << "Unknown graph type " << ((uint32_t)g_type);
  return nullptr;
}

template <template <typename, bool> typename OpType, typename... Args>
inline TraverseOperator* newTraverseOp(GraphType g_type, Args... args) {
  switch (g_type) {
  case GraphType::Normal:
    return new OpType<Graph, true>(std::forward<Args>(args)...);
  case GraphType::Partitioned:
    return new OpType<ReorderedPartitionedGraph, true>(std::forward<Args>(args)...);
  case GraphType::GraphView:
    return new OpType<GraphView<GraphPartitionBase>, true>(std::forward<Args>(args)...);
  case GraphType::BipartiteGraphView:
    return new OpType<GraphView<BipartiteGraph>, false>(std::forward<Args>(args)...);
  }
  LOG(FATAL) << "Unknown graph type " << ((uint32_t)g_type);
  return nullptr;
}

template <typename G>
inline constexpr bool sensitive_to_hint = !std::is_same<G, Graph>::value;

}  // namespace circinus

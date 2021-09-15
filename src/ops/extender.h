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

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "ops/operator.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

class Extender : public Operator {
 protected:
  const QueryGraph* query_graph_;
  const QueryVertexID query_vertex_;
  const std::vector<QueryVertexID> pivot_vertices_;
  std::string name_;
  VertexID filter_size_ = 0;
};

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

#include <string>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"

namespace circinus {

class QueryConfig {
 public:
  explicit QueryConfig(const std::string& config_str = "") {
    // FIXME(tatiana): parse string
  }
  std::string matching_order;
  std::string candidate_pruning_strategy = "cfl";
  std::string compression_strategy = "dynamic";
  bool use_auxiliary_index = false;
  bool use_partitioned_graph = true;
  std::string output = "count";
};

struct QueryContext {
  QueryGraph query_graph;
  QueryConfig query_config;
  Graph* data_graph;

  QueryContext(QueryGraph&& q, QueryConfig&& config, Graph* g)
      : query_graph(std::move(q)), query_config(std::move(config)), data_graph(g) {}

  void operator=(QueryContext&& context) {
    query_graph = std::move(context.query_graph);
    query_config = std::move(context.query_config);
    data_graph = context.data_graph;
  }
};

}  // namespace circinus

#pragma once

#include <string>
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

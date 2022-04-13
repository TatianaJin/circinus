#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/logical/filter/filter.h"

namespace circinus {

class NeighborhoodFilter;

class LogicalGQLFilter : public LogicalNeighborhoodFilter {
 public:
  explicit LogicalGQLFilter(const QueryGraph* query_graph);

  std::vector<std::unique_ptr<NeighborhoodFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                       ExecutionConfig& exec) override;
};

}  // namespace circinus

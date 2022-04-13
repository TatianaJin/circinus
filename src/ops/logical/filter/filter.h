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

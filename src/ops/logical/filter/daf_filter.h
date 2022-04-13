#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "ops/logical/filter/filter.h"

namespace circinus {

class NeighborhoodFilter;  // forward declaration

class LogicalDAFFilter : public LogicalNeighborhoodFilter {
 private:
  QueryVertexID start_vertex_;
  std::vector<TreeNode> tree_;
  std::vector<QueryVertexID> bfs_order_;

 public:
  QueryVertexID getStartVertex(const GraphMetadata& metadata, const QueryGraph* query_graph,
                               const std::vector<VertexID>& candidate_size);

  LogicalDAFFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                   const std::vector<VertexID>& candidate_size);

  const std::vector<QueryVertexID>& getBfsOrder() const { return bfs_order_; }

  std::vector<std::unique_ptr<NeighborhoodFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                       ExecutionConfig& exec) override;
};

}  // namespace circinus

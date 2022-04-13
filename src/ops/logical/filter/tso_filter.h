#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "ops/logical/filter/filter.h"

namespace circinus {

class NeighborhoodFilter;  // forward declaration

class LogicalTSOFilter : public LogicalNeighborhoodFilter {
 private:
  QueryVertexID start_vertex_;
  std::vector<TreeNode> tree_;
  std::vector<QueryVertexID> dfs_order_;

 public:
  LogicalTSOFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                   const std::vector<VertexID>& candidate_size);

  QueryVertexID getStartVertex(const GraphMetadata& metadata, const QueryGraph* query_graph_,
                               const std::vector<VertexID>& candidate_size);

  const std::vector<QueryVertexID>& getDfsOrder() const { return dfs_order_; }

  const std::vector<TreeNode>& getTree() const { return tree_; }

  std::vector<std::unique_ptr<NeighborhoodFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                       ExecutionConfig& exec) override;
};

}  // namespace circinus

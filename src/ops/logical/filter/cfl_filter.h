#pragma once

#include <memory>
#include <vector>

#include "algorithms/k_core.h"
#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/tree_node.h"
#include "ops/logical/filter/filter.h"

namespace circinus {

class NeighborhoodFilter;  // forward declaration
class Extender;            // forward declaration

class LogicalCFLFilter : public LogicalNeighborhoodFilter {
 private:
  QueryVertexID start_vertex_;
  std::vector<TreeNode> bfs_tree_;
  std::vector<QueryVertexID> bfs_order_;
  std::vector<QueryVertexID> level_offset_;
  uint32_t level_num_;
  TwoCoreSolver two_core_solver_;

 public:
  LogicalCFLFilter(const GraphMetadata& metadata, const QueryGraph* query_graph,
                   const std::vector<VertexID>& candidate_size, QueryVertexID seed_qv);

  /* Generating the BFS tree by Online query.
   */
  LogicalCFLFilter(const QueryGraph* query_graph, QueryVertexID seed_qv);

  void generateBFSTree(const QueryGraph* query_graph);

  std::vector<std::unique_ptr<NeighborhoodFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                       ExecutionConfig& exec) override;

  // std::vector<std::unique_ptr<Extender>> toPhysicalExtenders(const GraphMetadata& metadata,
  //                                                            ExecutionConfig& exec) override;

  QueryVertexID getStartVertex(const GraphMetadata& metadata, const QueryGraph* query_graph,
                               const std::vector<VertexID>& candidate_size);

  std::vector<QueryVertexID> getTopThree(const GraphMetadata& metadata, const QueryGraph* q);

  QueryVertexID getStartVertex(const std::vector<QueryVertexID>& query_vertices,
                               const std::vector<VertexID>& cardinality, const QueryGraph& q,
                               const GraphMetadata& metadata);

  const TwoCoreSolver& getTwoCoreSolver() const { return two_core_solver_; }

  const std::vector<TreeNode>& getTree() const { return bfs_tree_; }

  const std::vector<QueryVertexID>& getBfsOrder() const { return bfs_order_; }
};

}  // namespace circinus

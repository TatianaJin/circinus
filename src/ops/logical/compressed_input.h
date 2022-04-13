#pragma once

#include <memory>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

#include "glog/logging.h"

namespace circinus {

// forward declaration
class InputOperator;

class LogicalCompressedInputOperator {
 protected:
  const QueryVertexID input_query_vertex_;
  const bool inputs_are_keys_;

 public:
  LogicalCompressedInputOperator(const QueryVertexID input_query_vertex, const bool inputs_are_keys)
      : input_query_vertex_(input_query_vertex), inputs_are_keys_(inputs_are_keys) {}

  virtual ~LogicalCompressedInputOperator() {}

  virtual std::unique_ptr<InputOperator> toPhysicalOperators() const;
};

class PartitionedLogicalCompressedInputOperator : public LogicalCompressedInputOperator {
 protected:
  std::vector<std::pair<QueryVertexID, QueryVertexID>> qv_pivots_;

 public:
  PartitionedLogicalCompressedInputOperator(const QueryGraph* query_graph, bool inputs_are_keys,
                                            const std::vector<QueryVertexID>& matching_order,
                                            const std::vector<QueryVertexID>& partitioning_qvs);

  virtual ~PartitionedLogicalCompressedInputOperator() {}

  std::unique_ptr<InputOperator> toPhysicalOperators() const override;
};

}  // namespace circinus

#pragma once

#include <vector>

#include "graph/query_graph.h"

namespace circinus {

class TwoCoreSolver {
 private:
  const QueryGraph* graph_;
  std::vector<int> core_table_;

 public:
  inline static bool isInCore(const std::vector<int>& core_table, QueryVertexID v) { return core_table[v] > 1; }

  explicit TwoCoreSolver(const QueryGraph* graph);

  inline bool isInCore(QueryVertexID v) const { return core_table_[v] > 1; }

  inline uint32_t getCoreSize() {
    uint32_t count = 0;
    for (auto val : core_table_) {
      count += (val > 1);
    }
    return count;
  }

  const std::vector<int>& get2CoreTable() const { return core_table_; }

  std::vector<QueryVertexID> get2CoreVertices();

  inline QueryGraph extract2CoreSubgraph() { return graph_->getInducedSubgraph(get2CoreVertices()); }
};

}  // namespace circinus

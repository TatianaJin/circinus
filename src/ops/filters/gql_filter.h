#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters/filter.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

namespace circinus {

class GQLFilter : public NeighborhoodFilter {
 private:
  bool verify(const GraphBase* data_graph, const VertexID data_vertex,
              std::vector<std::vector<VertexID>>* candidates) const;

 public:
  GQLFilter(ExecutionConfig& conf, const QueryGraph* query_graph, QueryVertexID query_vertex)
      : NeighborhoodFilter(conf, query_graph, query_vertex) {}

  void filter(const GraphBase* data_graph, std::vector<std::vector<VertexID>>* candidates,
              FilterContext* ctx) const override;
};

}  // namespace circinus

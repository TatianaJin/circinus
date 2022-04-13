#pragma once

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/logical/filter/filter.h"
#include "utils/hashmap.h"

namespace circinus {

class LocalFilter;  // forward declaration

class LogicalNLFFilter : public LogicalLocalFilter {
 private:
  std::vector<unordered_map<LabelID, uint32_t>> out_neighbor_label_frequency_;
  std::vector<unordered_map<LabelID, uint32_t>> in_neighbor_label_frequency_;
  bool maximum_neighbor_degree_filter_ = false;
  std::vector<uint32_t> maximum_neighbor_degrees_;

 public:
  LogicalNLFFilter(const QueryGraph& q, const std::vector<QueryVertexID>& query_vertices,
                   bool maximum_neighbor_degree_filter)
      : out_neighbor_label_frequency_(query_vertices.size()),
        in_neighbor_label_frequency_(query_vertices.size()),
        maximum_neighbor_degree_filter_(maximum_neighbor_degree_filter) {
    for (uint32_t qi = 0; qi < query_vertices.size(); ++qi) {
      auto neighbors = q.getOutNeighbors(query_vertices[qi]);
      for (uint32_t i = 0; i < neighbors.second; ++i) {
        auto qv_label = q.getVertexLabel(neighbors.first[i]);
        if (qv_label == ALL_LABEL) continue;
        out_neighbor_label_frequency_[qi][qv_label] += 1;
      }
    }
    if (maximum_neighbor_degree_filter_) {
      maximum_neighbor_degrees_.resize(query_vertices.size(), 0);
      for (uint32_t qi = 0; qi < query_vertices.size(); ++qi) {
        auto neighbors = q.getOutNeighbors(query_vertices[qi]);
        for (uint32_t i = 0; i < neighbors.second; ++i) {
          maximum_neighbor_degrees_[qi] =
              std::max(maximum_neighbor_degrees_[qi], q.getVertexOutDegree(neighbors.first[i]));
        }
      }
    }
  }

  LogicalNLFFilter(const DirectedQueryGraph& q, const std::vector<QueryVertexID>& query_vertices,
                   bool maximum_neighbor_degree_filter)
      : out_neighbor_label_frequency_(query_vertices.size()), in_neighbor_label_frequency_(query_vertices.size()) {
    for (uint32_t qi = 0; qi < query_vertices.size(); ++qi) {
      auto neighbors = q.getOutNeighbors(query_vertices[qi]);
      for (uint32_t i = 0; i < neighbors.second; ++i) {
        auto qv_label = q.getVertexLabel(neighbors.first[i]);
        if (qv_label == ALL_LABEL) continue;
        out_neighbor_label_frequency_[qi][qv_label] += 1;
      }
      auto in_neighbors = q.getInNeighbors(query_vertices[qi]);
      for (uint32_t i = 0; i < in_neighbors.second; ++i) {
        auto qv_label = q.getVertexLabel(in_neighbors.first[i]);
        if (qv_label == ALL_LABEL) continue;
        in_neighbor_label_frequency_[qi][qv_label] += 1;
      }
    }
  }

  std::vector<std::unique_ptr<LocalFilter>> toPhysicalOperators(const GraphMetadata& metadata,
                                                                ExecutionConfig& exec) override;
};

}  // namespace circinus

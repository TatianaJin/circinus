#pragma once

#include <unordered_map>
#include <utility>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"
#include "ops/filters/local_filter.h"
#include "utils/hashmap.h"

namespace circinus {

class NLFFilter : public LocalFilter {
 private:
  unordered_map<LabelID, uint32_t> neighbor_label_frequency_;

 public:
  NLFFilter(const QueryGraph* query_graph, QueryVertexID query_vid);
  explicit NLFFilter(const unordered_map<LabelID, uint32_t>& neighbor_label_frequency)
      : neighbor_label_frequency_(neighbor_label_frequency) {
    DCHECK(!neighbor_label_frequency_.empty());
  }

  virtual ~NLFFilter() {}

  bool prune(const Graph& data_graph, VertexID v) const override;

  bool prune(const GraphPartition& g, VertexID v) const override;
};

class QuickNLFFilter : public LocalFilter {
 private:
  const unordered_map<LabelID, uint32_t> neighbor_label_frequency_;

 public:
  explicit QuickNLFFilter(const unordered_map<LabelID, uint32_t>& neighbor_label_frequency)
      : neighbor_label_frequency_(neighbor_label_frequency) {}

  bool prune(const Graph& data_graph, VertexID v) const override;

  bool prune(const GraphPartition& g, VertexID v) const override;
};

class MNDNLFFilter : public NLFFilter {
 private:
  const uint32_t maximum_neighbor_degree_;

 public:
  MNDNLFFilter(unordered_map<LabelID, uint32_t>&& neighbor_label_frequency, uint32_t maximum_neighbor_degree)
      : NLFFilter(std::move(neighbor_label_frequency)), maximum_neighbor_degree_(maximum_neighbor_degree) {}

  bool prune(const Graph& data_graph, VertexID v) const override;
};

}  // namespace circinus

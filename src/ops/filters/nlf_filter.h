// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

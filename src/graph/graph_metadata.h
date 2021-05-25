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

#include <vector>

#include "graph/graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

class GraphMetadata {
 private:
  bool directed_ = false;
  bool labeled_ = true;
  bool sorted_by_degree_ = false;
  bool sorted_by_label_ = false;
  bool has_label_index_ = true;
  VertexID n_vertices_ = 0;
  unordered_map<LabelID, VertexID> label_frequency_;
  // {degree, no. vertices with equal or greater out-degree}
  ordered_map<VertexID, VertexID> cardinality_by_out_degree_;
  // {degree, no. vertices with equal or greater in-degree}
  ordered_map<VertexID, VertexID> cardinality_by_in_degree_;

  /* for partitioned graph */
  VertexID vertex_offset_ = 0;
  std::vector<GraphMetadata> partitioned_metadata_;

 public:
  explicit GraphMetadata(const Graph& g) {
    collectBasicInfo(g);
    if (has_label_index_) {
      collectLabelFrequency(g);
    }
  }
  inline bool isDirected() const { return directed_; }
  inline bool isLabeled() const { return labeled_; }
  inline bool isSortedByDegree() const { return sorted_by_degree_; }
  inline bool isSortedByLabel() const { return sorted_by_label_; }
  inline bool hasLabelFrequency() const { return !label_frequency_.empty(); }
  inline bool hasOutDegreeFrequency() const { return !cardinality_by_out_degree_.empty(); }
  inline bool hasInDegreeFrequency() const { return !cardinality_by_in_degree_.empty(); }
  inline bool hasLabelIndex() const { return has_label_index_; }
  inline VertexID getNumVertices() const { return n_vertices_; }
  inline VertexID getLabelFrequency(LabelID label) const {
    auto pos = label_frequency_.find(label);
    if (pos == label_frequency_.end()) {
      return 0;
    }
    return pos->second;
  }

  inline VertexID getNumVerticesWithOutDegreeGE(VertexID degree) const {
    auto pos = cardinality_by_out_degree_.lower_bound(degree);
    if (pos == cardinality_by_out_degree_.end()) {
      return 0;
    }
    if (pos->first == degree) {
      return pos->second;
    }
    --pos;  // return an upper bound
    return pos->second;
  }

  inline VertexID getNumVerticesWithInDegreeGE(VertexID degree) const {
    auto pos = cardinality_by_in_degree_.lower_bound(degree);
    if (pos == cardinality_by_in_degree_.end()) {
      return 0;
    }
    if (pos->first == degree) {
      return pos->second;
    }
    --pos;  // return an upper bound
    return pos->second;
  }

  inline uint32_t numPartitions() const { return partitioned_metadata_.empty() ? 1 : partitioned_metadata_.size(); }
  const GraphMetadata& getPartition(uint32_t idx) const { return partitioned_metadata_[idx]; }

  // FIXME(tatiana)
  void collectBasicInfo(const Graph& graph) {
    n_vertices_ = graph.getNumVertices();
    labeled_ = graph.getNumLabels() > 1;
  }

  void collectDegreeFrequency(const Graph& graph) {}
  void collectLabelFrequency(const Graph& g) {
    for (auto label : g.getLabels()) {
      label_frequency_.insert({label, g.getVertexCardinalityByLabel(label)});
    }
  }
};

}  // namespace circinus

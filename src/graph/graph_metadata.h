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
#include <vector>

#include "graph/graph.h"
#include "graph/partitioned_graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

class GraphMetadata {
 private:
  bool directed_ = false;
  bool labeled_ = true;
  bool sorted_by_degree_ = false;  // this refers to secondary sorting
  bool sorted_by_label_ = false;   // this refers to primary sorting
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
  GraphMetadata(bool directed, bool labeled, bool sorted_by_degree, bool sorted_by_label, bool has_label_index,
                VertexID n_vertices, VertexID vertex_offset = 0)
      : directed_(directed),
        labeled_(labeled),
        sorted_by_degree_(sorted_by_degree),
        sorted_by_label_(sorted_by_label),
        has_label_index_(has_label_index),
        n_vertices_(n_vertices),
        vertex_offset_(vertex_offset) {}

  explicit GraphMetadata(const Graph& g, bool use_label_index = true)
      : GraphMetadata(false, g.getNumLabels() > 1, false, false, use_label_index, g.getNumVertices()) {
    if (use_label_index) {
      collectLabelFrequency(g);
    }
  }

  explicit GraphMetadata(const ReorderedPartitionedGraph& g, bool sorted_by_degree = true)
      : GraphMetadata(false, g.getNumLabels() > 1, false, g.getNumPartitions() == 1, false, g.getNumVertices()) {
    if (g.getNumPartitions() > 1) {
      partitioned_metadata_.reserve(g.getNumPartitions());
      for (uint32_t partition = 0; partition < g.getNumPartitions(); ++partition) {
        auto[start, end] = g.getPartitionRange(partition);
        partitioned_metadata_.emplace_back(directed_, labeled_, sorted_by_degree, true, false, end - start, start);
        // set label frequency in metadata
        auto& pm = partitioned_metadata_.back();
        auto& offsets = g.getLabelOffsetsInPartition(partition);
        for (auto& pair : offsets) {
          pm.label_frequency_.insert({pair.first, pair.second.second - pair.second.first});
          if (label_frequency_.find(pair.first) != label_frequency_.end()) {
            label_frequency_[pair.first] += pair.second.second - pair.second.first;
          } else {
            label_frequency_.insert({pair.first, pair.second.second - pair.second.first});
          }
        }
      }
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

  void collectDegreeFrequency(const GraphBase& graph) {
    // TODO(tatiana): see if needed
  }

  void collectLabelFrequency(const GraphBase& g) {
    for (auto label : g.getLabels()) {
      label_frequency_.insert({label, g.getVertexCardinalityByLabel(label)});
    }
  }
};

}  // namespace circinus

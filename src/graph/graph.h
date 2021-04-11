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

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"

namespace circinus {

class Graph {
 private:
  VertexID n_vertices_ = 0;
  EdgeID n_edges_ = 0;
  VertexID max_degree_ = 0;

  std::vector<EdgeID> vlist_;    // size n_vertices_ + 1, { i: the id of the first edge of vertex i }
  std::vector<VertexID> elist_;  // size n_edges_, { i : the destination vertex id of edge i}
  std::vector<LabelID> labels_;  // size n_vertices_, { i : the label of vertex i }

  std::unordered_map<LabelID, uint32_t> vertex_cardinality_by_label_;
  std::unordered_map<LabelID, std::vector<VertexID>> vertex_ids_by_label_;

 public:
  /**
   * @param path The file to load graph (undirected).
   * TODO(tatiana): support directed graph
   *
   * File format
   * t n_vertices n_edges
   * v label_id degree     => n_vertices lines, with continuous vertex id from 0
   * e src_id dst_id 0     => n_edges lines, put the end vertex with smaller id first for each edge
   */
  explicit Graph(const std::string& path);
  Graph() {}

  /// graph metadata
  inline VertexID getNumVertices() const { return n_vertices_; }
  inline EdgeID getNumEdges() const { return n_edges_; }
  inline VertexID getGraphMaxDegree() const { return max_degree_; }
  inline std::vector<LabelID> getLabels() const {
    std::vector<LabelID> labels;
    labels.reserve(vertex_ids_by_label_.size());
    for (auto& pair : vertex_ids_by_label_) {
      labels.push_back(pair.first);
    }
    return labels;
  }

  inline uint64_t getVertexCardinalityByLabel(LabelID label) const {
    auto pos = vertex_cardinality_by_label_.find(label);
    if (pos != vertex_cardinality_by_label_.end()) {
      return pos->second;
    }
    return 0;
  }

  /// vertex accessers
  inline VertexID getVertexOutDegree(VertexID id) const { return vlist_[id + 1] - vlist_[id]; }
  inline LabelID getVertexLabel(VertexID id) const { return labels_[id]; }

  inline const std::vector<VertexID>* getVerticesByLabel(LabelID lid) const {
    if (vertex_ids_by_label_.count(lid) == 0) {
      return nullptr;
    }
    return &vertex_ids_by_label_.at(lid);
  }

  inline const uint32_t getNumVerticesByLabel(LabelID lid) const {
    if (vertex_ids_by_label_.count(lid) == 0) {
      return 0;
    }
    return vertex_ids_by_label_.at(lid).size();
  }

  /** * @returns a pair { starting neighbor pointer, out degree } */
  inline std::pair<const VertexID*, uint32_t> getOutNeighbors(VertexID id) const {
    DCHECK_LT(id, vlist_.size() - 1);
    return std::make_pair(&elist_[vlist_[id]], vlist_[id + 1] - vlist_[id]);
  }

  /// persistence
  void DumpToFile(const std::string& path) const;
};

}  // namespace circinus

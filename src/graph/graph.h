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

#include "./metis.h"
#include "graph/graph_base.h"
#include "graph/types.h"
#include "graph/vertex_set_view.h"
#include "utils/hashmap.h"

namespace circinus {

class Graph : public GraphBase {
 protected:
  std::vector<LabelID> labels_;  // size n_vertices_, { i : the label of vertex i }
  unordered_map<LabelID, std::vector<VertexID>> vertex_ids_by_label_;

 public:
  /**
   * @param path The file to load graph (undirected).
   *
   * File format
   * t n_vertices n_edges
   * v label_id degree     => n_vertices lines, with continuous vertex id from 0
   * e src_id dst_id 0     => n_edges lines, put the end vertex with smaller id first for each edge
   */
  explicit Graph(const std::string& path, bool build_label_index = true);
  Graph() {}

  void buildLabelIndex();

  /* start of vertex-level accessers */

  /**
   * nbr_label Unused, added for the sake of GraphView-like interface.
   */
  inline VertexID getVertexOutDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx = 0) const {
    return getVertexOutDegree(id);
  }

  inline VertexID getVertexInDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx = 0) const {
    return getVertexOutDegree(id);
  }

  /* compatible interface for graph partition TODO(tatiana): unnecessary? */
  inline VertexID getVertexGlobalId(VertexID id) const { return id; }
  inline VertexID getVertexLocalId(VertexID id) const { return id; }

  inline LabelID getVertexLabel(VertexID id) const { return labels_[id]; }

  inline SingleRangeVertexSetView getOutNeighborsWithHint(VertexID id, LabelID nbr_label,
                                                          uint32_t graph_idx = 0) const {
    auto[first, size] = getOutNeighbors(id);
    return SingleRangeVertexSetView(first, size);
  }

  inline SingleRangeVertexSetView getInNeighborsWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx = 0) const {
    return getOutNeighborsWithHint(id, nbr_label, graph_idx);
  }

  /* end of vertex-level accessers */

  inline const std::vector<VertexID>* getVerticesByLabel(LabelID lid) const {
    DCHECK(!vertex_ids_by_label_.empty()) << "The label index is not constructed";
    if (vertex_ids_by_label_.count(lid) == 0) {
      return nullptr;
    }
    return &vertex_ids_by_label_.at(lid);
  }

  const auto& getVertexLabels() const { return labels_; }

  void clear() override {
    GraphBase::clear();
    labels_.clear();
    vertex_ids_by_label_.clear();
  }

  /// persistence
  void dumpToFile(const std::string& path) const override;

  std::pair<double, double> getMemoryUsage() const override;

 protected:
  void loadUndirectedGraphFromBinary(std::istream& input) override;
  void saveAsBinaryInner(std::ostream& output) const override;
};

class DirectedGraph : public Graph {
 public:
  VertexID getVertexInDegree(VertexID) const {
    // TODO(tatiana)
    return 0;
  }
};

}  // namespace circinus

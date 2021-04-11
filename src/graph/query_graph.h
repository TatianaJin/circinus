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
#include <utility>
#include <vector>
#include "glog/logging.h"

#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

struct QueryEdge {
  QueryVertexID src;
  QueryVertexID dst;

  QueryEdge(QueryVertexID _src, QueryVertexID _dst) : src(_src), dst(_dst) {}

  bool hasEndVertex(QueryVertexID v) const { return src == v || dst == v; }
};

/**
 * A QueryGraph is stored in the CSR format.
 * QueryGraph must have zero-indexed continuous vertex ids.
 */
class QueryGraph {
 private:
  QueryVertexID n_vertices_ = 0;
  EdgeID n_edges_ = 0;
  QueryVertexID max_degree_ = 0;

  std::vector<EdgeID> vlist_;         // size n_vertices_ + 1, { i: the id of the first edge of vertex i }
  std::vector<QueryVertexID> elist_;  // size n_edges_, { i : the destination vertex id of edge i}
  std::vector<LabelID> labels_;       // size n_vertices_, { i : the label of vertex i }
  unordered_map<LabelID, uint32_t> vertex_cardinality_by_label_;

 public:
  QueryGraph() {}

  /**
   * @param path The file to load query graph (undirected).
   * TODO(tatiana): support directed query graph
   *
   * File format
   * t n_vertices n_edges
   * v label_id degree     => n_vertices lines, with continuous vertex id from 0
   * e src_id dst_id       => n_edges lines, put the end vertex with smaller id first for each edge
   */
  explicit QueryGraph(const std::string& path);

  inline QueryVertexID getNumVertices() const { return n_vertices_; }
  inline QueryVertexID getGraphMaxDegree() const { return max_degree_; }
  inline uint32_t getVertexCardinalityByLabel(LabelID label) const {
    auto pos = vertex_cardinality_by_label_.find(label);
    return pos == vertex_cardinality_by_label_.end() ? 0 : pos->second;
  }
  inline QueryVertexID getVertexOutDegree(QueryVertexID id) const { return vlist_[id + 1] - vlist_[id]; }
  inline LabelID getVertexLabel(QueryVertexID id) const { return labels_[id]; }

  /** * @returns a pair { starting neighbor pointer, out degree } */
  inline std::pair<const QueryVertexID*, uint32_t> getOutNeighbors(QueryVertexID id) const {
    DCHECK_LT(id, vlist_.size() - 1);
    return std::make_pair(&elist_[vlist_[id]], vlist_[id + 1] - vlist_[id]);
  }

  /** extract induced subgraph on vertices, with the vertex index in vertices as the new vertex id */
  QueryGraph getInducedSubgraph(const std::vector<QueryVertexID>& vertices) const;
};

}  // namespace circinus

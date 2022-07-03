#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "./metis.h"
#include "glog/logging.h"

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
  using NeighborSet = SingleRangeVertexSetView;
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

  void loadUndirectedGraphEdgeList(const std::string& path);

  void buildLabelIndex();

  void reorderByDegree(bool ascending_degree = true);

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

  inline VertexID getVertexGlobalId(VertexID id) const { return id; }
  inline VertexID getVertexLocalId(VertexID id) const { return id; }

  inline LabelID getVertexLabel(VertexID id) const { return labels_[id]; }

  inline NeighborSet getOutNeighborsWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx = 0) const {
    auto[first, size] = getOutNeighbors(id);
    return NeighborSet(first, size);
  }

  inline NeighborSet getInNeighborsWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx = 0) const {
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
  VertexID getVertexInDegree(VertexID) const { return 0; }
};

}  // namespace circinus

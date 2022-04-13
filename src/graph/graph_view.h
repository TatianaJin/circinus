#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "graph/graph.h"
#include "graph/types.h"
#include "graph/vertex_set_view.h"

namespace circinus {

template <typename G>
class GraphView {
  /** Each graph corresponds to a query graph edge in `Expand` */
  std::vector<std::unique_ptr<G>> graphs_;  // owned, need to delete upon destruction

 public:
  using NeighborSet = std::conditional_t<std::is_same_v<Graph, G>, SingleRangeVertexSetView, VertexSetView>;
  explicit GraphView(std::vector<std::unique_ptr<G>>&& graphs) : graphs_(std::move(graphs)) {}

  /**
   * @param nbr_label Hint of the label of out neighbors. The results do not guarantee that the hint is used.
   */
  inline VertexID getVertexOutDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx) const {
    return graphs_[graph_idx]->getVertexOutDegreeWithHint(id, nbr_label);
  }

  inline VertexID getVertexInDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx) const {
    return graphs_[graph_idx]->getVertexInDegreeWithHint(id, nbr_label);
  }

  /**
   * @param nbr_label Hint of the label of out neighbors. The results do not guarantee that the hint is used.
   */
  inline NeighborSet getOutNeighborsWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx) const {
    return graphs_[graph_idx]->getOutNeighborsWithHint(id, nbr_label);
  }

  inline NeighborSet getInNeighborsWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx) const {
    return graphs_[graph_idx]->getInNeighborsWithHint(id, nbr_label);
  }
};

}  // namespace circinus

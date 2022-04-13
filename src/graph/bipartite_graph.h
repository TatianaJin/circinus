#pragma once

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/candidate_set_view.h"
#include "graph/graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

/** The projection of the data graph on two sets of vertices */
class BipartiteGraph
    : public Graph {  // Overwrite variable:vlist_,elist_,n_vertices_  function:getVertexOutDegree,getOutNeighbors
 private:
  unordered_map<VertexID, uint32_t> offset_by_vertex_;
  bool populated_ = false;
  QueryVertexID source_id_, destination_id_;
  uint64_t bipartite_graph_intersection_input_size_ = 0;
  uint64_t bipartite_graph_intersection_output_size_ = 0;

 public:
  using NeighborSet = VertexSetView;
  BipartiteGraph(VertexID id1, VertexID id2) : Graph(), source_id_(id1), destination_id_(id2) {}

  inline QueryVertexID getSourceId() const { return source_id_; }

  inline void populateGraph(const GraphBase* g, const std::vector<CandidateSetView>& candidate_sets) {
    DCHECK_LT(source_id_, candidate_sets.size());
    DCHECK_LT(destination_id_, candidate_sets.size());
    populateGraph(g, candidate_sets[source_id_], candidate_sets[destination_id_]);
  }

  /** Only populates the edges from vertices in candidate_set1 to vertices in candidate_set2 */
  void populateGraph(const GraphBase* g, const CandidateSetView& candidate_set1,
                     const CandidateSetView& candidate_set2) {
    if (populated_) return;

    populated_ = true;
    unordered_set<VertexID> vset(candidate_set2.begin(), candidate_set2.end());
    n_vertices_ = candidate_set1.size();
    vlist_.reserve(n_vertices_ + 1);
    vlist_.emplace_back(0);

    size_t i = 0;
    for (auto v1_id : candidate_set1) {
      offset_by_vertex_.insert({v1_id, i++});
      auto[dest_nodes, cnt] = g->getOutNeighbors(v1_id);

      for (uint32_t j = 0; j < cnt; ++j) {
        if (vset.find(dest_nodes[j]) != vset.end()) elist_.emplace_back(dest_nodes[j]);
      }
      vlist_.emplace_back(elist_.size());
    }
    bipartite_graph_intersection_input_size_ += candidate_set1.size() + candidate_set2.size();
    bipartite_graph_intersection_output_size_ += elist_.size();
    n_edges_ = elist_.size();
  }

  std::pair<uint64_t, uint64_t> getProfilePair() const {
    return {bipartite_graph_intersection_input_size_, bipartite_graph_intersection_output_size_};
  }

  inline uint32_t getOffset(VertexID id) const {
    CHECK_EQ(offset_by_vertex_.count(id), 1) << " id " << id;
    return offset_by_vertex_.at(id);
  }

  /**
   * @returns The original ids of neighbor vertices of the given vertex
   */
  inline NeighborSet getOutNeighborsWithHint(VertexID id, LabelID nbr_label = ALL_LABEL) const {
    DCHECK_EQ(offset_by_vertex_.count(id), 1);
    const uint32_t offset = offset_by_vertex_.at(id);
    return Graph::getOutNeighborsWithHint(offset, nbr_label);
  }

  inline VertexID getVertexOutDegree(VertexID id, LabelID nbr_label = ALL_LABEL) const {
    DCHECK_EQ(offset_by_vertex_.count(id), 1);
    const uint32_t offset = offset_by_vertex_.at(id);
    return Graph::getVertexOutDegree(offset);
  }
};

}  // namespace circinus

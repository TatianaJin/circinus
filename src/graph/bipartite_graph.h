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
#include <unordered_set>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

/** The projection of the data graph on two sets of vertices */
class BipartiteGraph : public Graph {  // only use variable:vlist_,elist_  function:getVertexOutDegree,getOutNeighbors
 private:
  unordered_map<VertexID, uint32_t> offset_by_vertex_;
  bool populated_ = 0;
  QueryVertexID source_id_, destination_id_;
  uint64_t bipartite_graph_intersection_input_size_ = 0;
  uint64_t bipartite_graph_intersection_output_size_ = 0;

 public:
  BipartiteGraph(VertexID id1, VertexID id2) : Graph(), source_id_(id1), destination_id_(id2) {}

  void populateGraph(const Graph* g, const std::vector<std::vector<VertexID>>* candidate_sets) {
    populateGraph(g, (*candidate_sets)[source_id_], (*candidate_sets)[destination_id_]);
  }

  /** Only populates the edges from vertices in candidate_set1 to vertices in candidate_set2 */
  void populateGraph(const Graph* g, const std::vector<VertexID>& candidate_set1,
                     const std::vector<VertexID>& candidate_set2) {
    if (populated_) return;
    populated_ = 1;
    unordered_set<VertexID> vset(candidate_set2.begin(), candidate_set2.end());
    vlist_.emplace_back(0);
    for (size_t i = 0; i < candidate_set1.size(); ++i) {
      VertexID v1Id = candidate_set1[i];
      offset_by_vertex_.insert({v1Id, i});
      auto[dest_nodes, cnt] = g->getOutNeighbors(v1Id);
      for (uint32_t j = 0; j < cnt; ++j)
        if (vset.find(dest_nodes[j]) != vset.end()) elist_.emplace_back(dest_nodes[j]);
      vlist_.emplace_back(elist_.size());
    }
    bipartite_graph_intersection_input_size_ += candidate_set1.size() + candidate_set2.size();
    bipartite_graph_intersection_output_size_ += elist_.size();
  }

  std::pair<uint64_t, uint64_t> getProfilePair() const {
    return {bipartite_graph_intersection_input_size_, bipartite_graph_intersection_output_size_};
  }

  /**
   * @returns The original ids of neighbor vertices of the given vertex
   */
  inline std::pair<const VertexID*, uint32_t> getOutNeighbors(VertexID id, LabelID dummy = 0) const {
    DCHECK_EQ(offset_by_vertex_.count(id), 1);
    const uint32_t offset = offset_by_vertex_.at(id);
    return Graph::getOutNeighbors(offset);
  }

  inline VertexID getVertexOutDegree(VertexID id, LabelID dummy = 0) const {
    DCHECK_EQ(offset_by_vertex_.count(id), 1);
    const uint32_t offset = offset_by_vertex_.at(id);
    return Graph::getVertexOutDegree(offset);
  }
};

}  // namespace circinus

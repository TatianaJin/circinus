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

#include <utility>
#include <vector>

#include "graph/graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

class BipartiteGraph : public Graph {  // only use variable:vlist_,elist_  function:getVertexOutDegree,getOutNeighbors
 private:
  unordered_map<VertexID, uint32_t> offset_by_vertex_;
  bool populated = 0;
  VertexID sourceId, destinationId;

 public:
  BipartiteGraph(VertexID id1, VertexID id2) : Graph(), sourceId(id1), destinationId(id2) {}

  void populateGraph(const Graph* g, std::vector<std::vector<VertexID>>* candidate_sets) {
    populateGraph(g, (*candidate_sets)[sourceId], (*candidate_sets)[destinationId]);
  }

  void populateGraph(const Graph* g, std::vector<VertexID>& candidate_set1, std::vector<VertexID>& candidate_set2) {
    if (populated) return;
    populated = 1;
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
  }

  inline std::pair<const VertexID*, uint32_t> getOutNeighbors(
      VertexID id) const {  // populated must be 1 now, but dont add if-else here for performance
    const uint32_t offset = offset_by_vertex_.at(id);
    return std::make_pair(&elist_[vlist_[offset]], vlist_[offset + 1] - vlist_[offset]);
  }
  inline VertexID getVertexOutDegree(
      VertexID id) const {  // populated must be 1 now, but dont add if-else here for performance
    const uint32_t offset = offset_by_vertex_.at(id);
    return vlist_[offset + 1] - vlist_[offset];
  }
};

}  // namespace circinus

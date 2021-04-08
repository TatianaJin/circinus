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

 public:
  explicit BipartiteGraph(Graph g, std::vector<VertexID> candidate_set1, std::vector<VertexID> candidate_set2)
      : Graph() {
    unordered_set<VertexID> vset;
    vlist_.emplace_back(0);
    for (size_t i = 0; i < candidate_set1.size(); ++i) {
      VertexID v1Id = candidate_set1[i];
      offset_by_vertex_.insert({v1Id, i});
      auto[dest_nodes, cnt] = g.getOutNeighbors(v1Id);
      vset.clear();
      for (uint32_t j = 0; j < cnt; ++j) vset.insert(dest_nodes[j]);
      for (VertexID v2Id : candidate_set2)
        if (vset.find(v2Id) != vset.end()) elist_.emplace_back(v2Id);
      vlist_.emplace_back(elist_.size());
    }
  }
  inline std::pair<const VertexID*, uint32_t> getOutNeighbors(VertexID id) const {
    const uint32_t offset = offset_by_vertex_.at(id);
    return std::make_pair(&elist_[vlist_[offset]], vlist_[offset + 1] - vlist_[offset]);
  }
  inline VertexID getVertexOutDegree(VertexID id) const {
    const uint32_t offset = offset_by_vertex_.at(id);
    return vlist_[offset + 1] - vlist_[offset];
  }
};

}  // namespace circinus

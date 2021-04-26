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
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/query_graph.h"
#include "ops/traverse_operator.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandVertexOperator : public TraverseOperator {
 protected:
  std::vector<QueryVertexID> parents_;
  QueryVertexID target_vertex_;
  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;

 public:
  ExpandVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                       const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                       const std::vector<uint32_t>& same_label_key_indices,
                       const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                       SubgraphFilter* subgraph_filter)
      : TraverseOperator(same_label_key_indices, same_label_set_indices, set_pruning_threshold, subgraph_filter),
        parents_(parents),
        target_vertex_(target_vertex),
        query_vertex_indices_(query_vertex_indices) {}

  virtual ~ExpandVertexOperator() {}

  const auto& getQueryVertexIndices() const { return query_vertex_indices_; }

  std::vector<BipartiteGraph*> computeBipartiteGraphs(
      const Graph* g, const std::vector<std::vector<VertexID>>& candidate_sets) override {
    std::vector<BipartiteGraph*> ret;
    ret.reserve(parents_.size());
    for (auto parent_vertex : parents_) {
      ret.emplace_back(new BipartiteGraph(parent_vertex, target_vertex_));
      ret.back()->populateGraph(g, &candidate_sets);
    }
    return ret;
  }

 protected:
  inline void toStringInner(std::stringstream& ss) const {
    for (auto parent : parents_) {
      DCHECK_EQ(query_vertex_indices_.count(parent), 1);
      ss << ' ' << parent;
    }
    DCHECK_EQ(query_vertex_indices_.count(target_vertex_), 1);
    ss << " -> " << target_vertex_;
    if (candidates_ != nullptr) ss << " (" << candidates_->size() << ")";
  }
};

}  // namespace circinus

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

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "algorithms/intersect.h"
#include "graph/compressed_subgraphs.h"
#include "graph/query_graph.h"
#include "ops/expand_vertex_traverse_context.h"
#include "ops/traverse_operator.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandVertexOperator : public TraverseOperator {
 protected:
  std::vector<QueryVertexID> parents_;
  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;

 public:
  ExpandVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                       const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                       const std::vector<uint32_t>& same_label_key_indices,
                       const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                       SubgraphFilter* subgraph_filter)
      : TraverseOperator(target_vertex, same_label_key_indices, same_label_set_indices, set_pruning_threshold,
                         subgraph_filter),
        parents_(parents),
        query_vertex_indices_(query_vertex_indices) {}

  virtual ~ExpandVertexOperator() {}

  const auto& getQueryVertexIndices() const { return query_vertex_indices_; }

  std::unique_ptr<TraverseContext> initTraverseContext(const std::vector<CompressedSubgraphs>* inputs,
                                                       const void* graph, uint32_t start, uint32_t end,
                                                       QueryType profile) const override {
    auto ret = std::make_unique<ExpandVertexTraverseContext>(inputs, graph, start, end, parents_.size());
    ret->query_type = profile;
    return ret;
  }

  std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<CandidateSetView>& candidate_sets) override {
    std::vector<std::unique_ptr<BipartiteGraph>> ret;
    ret.reserve(parents_.size());
    for (auto parent_vertex : parents_) {
      ret.emplace_back(std::make_unique<BipartiteGraph>(parent_vertex, target_vertex_));
      ret.back()->populateGraph(g, candidate_sets);
    }
    return ret;
  }

  std::vector<std::unique_ptr<GraphPartitionBase>> computeGraphPartitions(
      const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const override {
    std::vector<std::unique_ptr<GraphPartitionBase>> ret;
    ret.reserve(parents_.size());
    for (auto parent_vertex : parents_) {
      ret.emplace_back(GraphPartitionBase::createGraphPartition(candidate_scopes[parent_vertex],
                                                                candidate_scopes[target_vertex_], g));
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

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

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/intersect.h"
#include "graph/compressed_subgraphs.h"
#include "graph/query_graph.h"
#include "ops/expand_vertex_traverse_context.h"
#include "ops/traverse_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandVertexOperator : public TraverseOperator {
 protected:
  std::vector<QueryVertexID> parents_;
  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;
  std::vector<uint32_t> parent_indices_;

 public:
  ExpandVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                       const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                       const std::vector<uint32_t>& same_label_key_indices,
                       const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                       std::unique_ptr<SubgraphFilter>&& subgraph_filter)
      : TraverseOperator(target_vertex, same_label_key_indices, same_label_set_indices, set_pruning_threshold,
                         std::move(subgraph_filter)),
        parents_(parents),
        query_vertex_indices_(query_vertex_indices),
        parent_indices_(parents.size()) {
    for (uint32_t i = 0; i < parents_.size(); ++i) {
      parent_indices_[i] = query_vertex_indices_.at(parents_[i]);
    }
  }

  virtual ~ExpandVertexOperator() {}

  bool extend_vertex() const override { return true; }

  const auto& getQueryVertexIndices() const { return query_vertex_indices_; }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    return std::make_unique<ExpandVertexTraverseContext>(candidates, graph, outputs, profile, parents_.size());
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

  template <typename G, QueryType profile, bool intersect_candidates>
  void expandFromParents(const CompressedSubgraphs& input, const G* data_graph, ExpandVertexTraverseContext* ctx,
                         const std::vector<uint32_t>& parent_indices, const unordered_set<VertexID>& exceptions,
                         std::vector<VertexID>* targets) const {
    using NeighborSet = typename G::NeighborSet;
    std::vector<NeighborSet> neighbor_sets;
    neighbor_sets.reserve(parent_indices.size());
    for (uint32_t i = 0; i < parent_indices.size(); ++i) {
      DCHECK_LT(parent_indices[i], input.getNumKeys());
      uint32_t key_vid = input.getKeyVal(parent_indices[i]);
      neighbor_sets.push_back(data_graph->getOutNeighborsWithHint(key_vid, target_label_, i));
    }
    std::sort(neighbor_sets.begin(), neighbor_sets.end(),
              [](const NeighborSet& a, const NeighborSet& b) { return a.size() < b.size(); });
    if (neighbor_sets.front().size() == 0) return;

    bool intersect_candidates_first = intersect_candidates && (!isProfileCandidateSIEffect(profile)) &&
                                      ctx->getCandidateSet()->size() <= neighbor_sets.front().size();
    uint32_t loop_idx = 2;
    uint32_t si_input_size = neighbor_sets.front().size();
    if (intersect_candidates_first) {
      loop_idx = 1;
      intersect(neighbor_sets.front(), *ctx->getCandidateSet(), targets, exceptions);
      si_input_size += ctx->getCandidateSet()->size();
    } else {
      intersect(neighbor_sets[0], neighbor_sets[1], targets, exceptions);
      si_input_size += neighbor_sets[1].size();
    }
    if
      constexpr(isProfileMode(profile)) { ctx->updateIntersectInfo(si_input_size, targets->size()); }
    if (targets->empty()) {
      return;
    }
    {  // enforce partial order
      filterTargets(targets, input);
      if (targets->empty()) return;
    }
    for (uint32_t i = loop_idx; i < parent_indices.size(); ++i) {
      auto si_input_size = targets->size() + neighbor_sets[i].size();
      (void)si_input_size;
      intersectInplace(*targets, neighbor_sets[i], targets);
      if
        constexpr(isProfileMode(profile)) { ctx->updateIntersectInfo(si_input_size, targets->size()); }
      if (targets->empty()) {
        return;
      }
    }
    if (intersect_candidates && !intersect_candidates_first) {
      auto new_key_size = targets->size();
      (void)new_key_size;
      intersectInplace(*targets, *ctx->getCandidateSet(), targets);
      if
        constexpr(isProfileMode(profile)) {
          ctx->updateIntersectInfo(new_key_size + ctx->getCandidateSet()->size(), targets->size());
          ctx->candidate_si_diff += new_key_size - targets->size();
        }
    } else if (!intersect_candidates) {
      degreeFilter(targets, target_degree_, data_graph);
    }
  }

 protected:
  inline void toStringInner(std::stringstream& ss) const {
    for (auto parent : parents_) {
      DCHECK_EQ(query_vertex_indices_.count(parent), 1);
      ss << ' ' << parent;
    }
    DCHECK_EQ(query_vertex_indices_.count(target_vertex_), 1);
    ss << " -> " << target_vertex_ << ':' << target_label_;
  }
};

}  // namespace circinus

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

#include <chrono>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/bipartite_graph.h"
#include "graph/candidate_set_view.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/graph_partition.h"
#include "graph/partitioned_graph.h"
#include "graph/types.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/filters/target_filter.h"
#include "ops/operator.h"
#include "ops/traverse_context.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class TraverseOperator : public Operator {
 protected:
  const QueryVertexID target_vertex_;
  LabelID target_label_ = ALL_LABEL;

  /* for non-repeated-vertex check and/or automorphism elimination */
  uint64_t set_pruning_threshold_ = ~0u;
  std::vector<uint32_t> same_label_key_indices_;
  std::vector<uint32_t> same_label_set_indices_;
  std::unique_ptr<SubgraphFilter> subgraph_filter_ = nullptr;
  std::unique_ptr<TargetFilter> target_filter_ = nullptr;

  std::vector<std::pair<bool, uint32_t>> matching_order_indices_;

 public:
  TraverseOperator(QueryVertexID target_vertex, std::unique_ptr<SubgraphFilter>&& filter)
      : target_vertex_(target_vertex), subgraph_filter_(std::move(filter)) {}

  TraverseOperator(QueryVertexID target_vertex, const std::vector<uint32_t>& same_label_key_indices,
                   const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                   std::unique_ptr<SubgraphFilter>&& filter = nullptr)
      : target_vertex_(target_vertex),
        set_pruning_threshold_(set_pruning_threshold),
        same_label_key_indices_(same_label_key_indices),
        same_label_set_indices_(same_label_set_indices),
        subgraph_filter_(std::move(filter)) {}
  virtual ~TraverseOperator() {}

  /* setters */
  inline void setMatchingOrderIndices(std::vector<std::pair<bool, uint32_t>>&& matching_order_indices) {
    matching_order_indices_ = std::move(matching_order_indices);
  }
  inline void setTargetLabel(LabelID l) { target_label_ = l; }

  inline void filterTargets(std::vector<VertexID>* targets, const CompressedSubgraphs& group) const {
    if (target_filter_ == nullptr) return;
    target_filter_->filter(targets, group);
  }

  /** @returns True if target is to be pruned. */
  inline bool filterTarget(VertexID target, const CompressedSubgraphs& group) const {
    if (target_filter_ == nullptr) return false;
    return target_filter_->filter(target, group);
  }

  // TODO(tatiana): add degree filter for non-candidate-intersection cases
  void addFilters(const std::vector<std::pair<bool, uint32_t>>& lt_constraints,
                  const std::vector<std::pair<bool, uint32_t>>& gt_constraints) {
    std::vector<uint32_t> expand_lt_conditions, expand_gt_conditions;      // key indices
    std::vector<uint32_t> subgraph_lt_conditions, subgraph_gt_conditions;  // set indices
    for (auto& p : lt_constraints) {
      if (p.first) {
        expand_lt_conditions.push_back(p.second);
      } else {
        subgraph_lt_conditions.push_back(p.second);
      }
    }
    for (auto& p : gt_constraints) {
      if (p.first) {
        expand_gt_conditions.push_back(p.second);
      } else {
        subgraph_gt_conditions.push_back(p.second);
      }
    }
    if (!subgraph_lt_conditions.empty() || !subgraph_gt_conditions.empty()) {
      if (subgraph_filter_ == nullptr) {
        subgraph_filter_ = SubgraphFilter::newPartialOrderSubgraphFilter(
            std::move(subgraph_lt_conditions), std::move(subgraph_gt_conditions), matching_order_indices_.back());
      } else {
        subgraph_filter_ = subgraph_filter_->addPartialOrderSubgraphFilter(
            std::move(subgraph_lt_conditions), std::move(subgraph_gt_conditions), matching_order_indices_.back());
      }
    }
    if (expand_lt_conditions.empty() && expand_gt_conditions.empty()) return;
    target_filter_ = TargetFilter::newTargetFilter(std::move(expand_lt_conditions), std::move(expand_gt_conditions));
  }

  /* getters */
  inline const auto& getSameLabelKeyIndices() const { return same_label_key_indices_; }
  inline const auto& getSameLabelSetIndices() const { return same_label_set_indices_; }
  inline auto getSetPruningThreshold() const { return set_pruning_threshold_; }
  inline QueryVertexID getTargetQueryVertex() const { return target_vertex_; }
  inline const auto& getMatchingOrderIndices() const { return matching_order_indices_; }
  inline auto getTargetLabel() const { return target_label_; }
  TargetFilter* getTargetFilter() const { return target_filter_.get(); }

  inline bool filter(CompressedSubgraphs& subgraphs) const {
    DCHECK(subgraph_filter_ != nullptr);
    return subgraph_filter_->filter(subgraphs);
  }

  inline uint32_t filter(std::vector<CompressedSubgraphs>& subgraphs, uint32_t start, uint32_t end) const {
    DCHECK(subgraph_filter_ != nullptr);
    return subgraph_filter_->filter(subgraphs, start, end);
  }

  // for backward compabilitity of tests
  inline std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<std::vector<VertexID>>& candidate_sets) {
    std::vector<CandidateSetView> views(candidate_sets.begin(), candidate_sets.end());
    return computeBipartiteGraphs(g, views);
  }

  virtual std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<CandidateSetView>& candidate_sets) = 0;

  virtual std::vector<std::unique_ptr<GraphPartitionBase>> computeGraphPartitions(
      const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const = 0;

  virtual std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap = nullptr) const = 0;

  virtual std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const = 0;

  virtual uint32_t expand(uint32_t cap, TraverseContext* ctx) const = 0;

  virtual bool extend_vertex() const = 0;

  uint32_t expandAndProfile(uint32_t cap, TraverseContext* ctx) const {
    auto start = std::chrono::high_resolution_clock::now();
    auto n = expandAndProfileInner(cap, ctx);
    auto stop = std::chrono::high_resolution_clock::now();
    ctx->total_time_in_milliseconds += toMilliseconds(start, stop);
    if (false) {
      uint32_t size = 0;
      auto offset = ctx->getOutputSize() - n;
      auto& outputs = *(ctx->getOutputs());
      for (uint32_t i = 0; i < n; ++i) {
        if ((*(ctx->getOutputs()))[offset + i].getNumIsomorphicSubgraphs(1) == 0) {
          // CHECK(false) << toString() << "\n\t" << (*outputs)[offset + i].toString();
          continue;
        }
        if (size != i) {
          outputs[offset + size] = std::move(outputs[offset + i]);
        }
        ++size;
      }
      ctx->popOutputs(n - size);
      n = size;
      std::ofstream ss("output_" + std::to_string(outputs.front().getNumVertices()), std::ofstream::app);
      for (auto& output : outputs) {
        output.logEnumerated(ss, matching_order_indices_);
      }
    }
    ctx->total_num_output_subgraphs +=
        getNumSubgraphs(*(ctx->getOutputs()), ctx->getOutputSize() - n, ctx->getOutputSize());
    ctx->total_output_size += n;
    return n;
  }

  std::string toProfileString(const ProfileInfo& info) const override {
    std::stringstream ss;
    ss << toString() << ',' << info.total_time_in_milliseconds << ',' << info.total_input_size << ','
       << info.total_output_size << ',' << info.total_num_input_subgraphs << ',' << info.total_num_output_subgraphs
       << ',' << info.intersection_count << ',' << info.total_intersection_input_size << ','
       << info.total_intersection_output_size << ',' << info.distinct_intersection_count << ','
       << info.candidate_si_diff;
#ifdef INTERSECTION_CACHE
    ss << ',' << info.cache_hit;
#endif
    return ss.str();
  }

 protected:
  virtual uint32_t expandAndProfileInner(uint32_t cap, TraverseContext* ctx) const = 0;
};

}  // namespace circinus

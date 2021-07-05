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
#include "ops/operator.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class TraverseContext : public ProfileInfo {
 private:
  uint32_t input_index_ = 0;
  uint32_t input_start_index_ = 0;
  uint32_t input_end_index_ = 0;
  const std::vector<CompressedSubgraphs>* current_inputs_ = nullptr;

 public:
  const void* current_data_graph = nullptr;
  std::vector<CompressedSubgraphs>* outputs;
  QueryType query_type;

  TraverseContext() {}
  TraverseContext(const std::vector<CompressedSubgraphs>* inputs, const void* data_graph)
      : TraverseContext(inputs, data_graph, 0, inputs->size()) {}

  TraverseContext(const std::vector<CompressedSubgraphs>* inputs, const void* data_graph, uint32_t input_index,
                  uint32_t input_end_index)
      : input_index_(input_index),
        input_start_index_(input_index),
        input_end_index_(input_end_index),
        current_inputs_(inputs),
        current_data_graph(data_graph) {}

  virtual ~TraverseContext() {}

  inline const CompressedSubgraphs& getCurrentInput() const { return (*current_inputs_)[input_index_]; }
  inline const CompressedSubgraphs& getPreviousInput() const { return (*current_inputs_)[input_index_ - 1]; }
  inline bool hasNextInput() const { return input_index_ < input_end_index_; }
  inline uint32_t getInputIndex() const { return input_index_; }
  inline uint32_t getTotalInputSize() const { return input_index_ - input_start_index_; }
  inline const auto getOutputs() const { return outputs; }

  inline void nextInput() { ++input_index_; }
  inline void clearOutputs() { outputs->clear(); }
};

class TraverseOperator : public Operator {
 protected:
  const CandidateSetView* candidates_ = nullptr;
  const QueryVertexID target_vertex_;

  /* for non-repeated-vertex check */
  uint64_t set_pruning_threshold_ = ~0u;
  std::vector<uint32_t> same_label_key_indices_;
  std::vector<uint32_t> same_label_set_indices_;
  SubgraphFilter* const subgraph_filter_ = nullptr;  // owned by the execution plan

  std::vector<std::pair<bool, uint32_t>> matching_order_indices_;

 public:
  TraverseOperator(QueryVertexID target_vertex, SubgraphFilter* filter)
      : target_vertex_(target_vertex), subgraph_filter_(filter) {}
  TraverseOperator(QueryVertexID target_vertex, const std::vector<uint32_t>& same_label_key_indices,
                   const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                   SubgraphFilter* filter = nullptr)
      : target_vertex_(target_vertex),
        set_pruning_threshold_(set_pruning_threshold),
        same_label_key_indices_(same_label_key_indices),
        same_label_set_indices_(same_label_set_indices),
        subgraph_filter_(filter) {}
  virtual ~TraverseOperator() {}

  /* setters */
  inline virtual void setCandidateSets(const CandidateSetView* candidates) { candidates_ = candidates; }
  inline void setMatchingOrderIndices(std::vector<std::pair<bool, uint32_t>>&& matching_order_indices) {
    matching_order_indices_ = std::move(matching_order_indices);
  }

  /* getters */
  inline const CandidateSetView* getCandidateSet() const { return candidates_; }
  inline const auto& getSameLabelKeyIndices() const { return same_label_key_indices_; }
  inline const auto& getSameLabelSetIndices() const { return same_label_set_indices_; }
  inline auto getSetPruningThreshold() const { return set_pruning_threshold_; }
  inline QueryVertexID getTargetQueryVertex() const { return target_vertex_; }
  inline const auto& getMatchingOrderIndices() const { return matching_order_indices_; }

  inline bool filter(const CompressedSubgraphs& subgraphs) const {
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

  virtual std::unique_ptr<TraverseContext> initTraverseContext(const std::vector<CompressedSubgraphs>* inputs,
                                                               const void* graph, uint32_t input_start,
                                                               uint32_t input_end) const {
    return std::make_unique<TraverseContext>(inputs, graph, input_start, input_end);
  }

  virtual uint32_t expand(uint32_t cap, TraverseContext* ctx) const = 0;

  uint32_t expandAndProfile(uint32_t cap, TraverseContext* ctx) const {
    auto start = std::chrono::high_resolution_clock::now();
    auto n = expandAndProfileInner(cap, ctx);
    auto stop = std::chrono::high_resolution_clock::now();
    ctx->total_time_in_milliseconds +=
        (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000000.0);
    {
      uint32_t size = 0;
      auto offset = ctx->outputs->size() - n;
      for (uint32_t i = 0; i < n; ++i) {
        if ((*(ctx->outputs))[offset + i].getNumIsomorphicSubgraphs(1) == 0) {
          // CHECK(false) << toString() << "\n\t" << (*outputs)[offset + i].toString();
          continue;
        }
        if (size != i) {
          (*(ctx->outputs))[offset + size] = std::move((*(ctx->outputs))[offset + i]);
        }
        ++size;
      }
      ctx->outputs->erase(ctx->outputs->begin() + (offset + size), ctx->outputs->end());
      n = size;
    }
    if (false) {
      std::ofstream ss("output_" + std::to_string(ctx->outputs->front().getNumVertices()), std::ofstream::app);
      for (auto& output : *(ctx->outputs)) {
        output.logEnumerated(ss, matching_order_indices_);
      }
    }
    ctx->total_num_output_subgraphs += getNumSubgraphs(*(ctx->outputs), ctx->outputs->size() - n, ctx->outputs->size());
    ctx->total_output_size += n;
    return n;
  }

  // FIXME(tatiana): profile info interface
  std::string toProfileString(TraverseContext* ctx) const {
    std::stringstream ss;
    ss << toString() << ',' << ctx->total_time_in_milliseconds << ',' << ctx->getTotalInputSize() << ','
       << ctx->total_output_size << ',' << ctx->total_num_input_subgraphs << ',' << ctx->total_num_output_subgraphs
       << ',' << ctx->intersection_count << ',' << ctx->total_intersection_input_size << ','
       << ctx->total_intersection_output_size << ',' << ctx->distinct_intersection_count;
    return ss.str();
  }

 protected:
  virtual uint32_t expandAndProfileInner(uint32_t cap, TraverseContext* ctx) const = 0;
};

}  // namespace circinus

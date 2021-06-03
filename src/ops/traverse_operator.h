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
#include "ops/operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

namespace circinus {

inline void removeExceptions(const VertexSetView& set, std::vector<VertexID>* res,
                             const unordered_set<VertexID>& except = {}) {
  for (auto vid : set) {
    if (except.count(vid) == 0) res->emplace_back(vid);
  }
}

struct TraverseContext {
  /* transient variables for recording the current inputs */
  const void* current_data_graph = nullptr;

  uint32_t current_inputs_size = 0;
  std::vector<CompressedSubgraphs>* outputs;

  /* for profiling */
  uint64_t total_input_size = 0;
  uint64_t total_output_size = 0;
  uint64_t total_num_input_subgraphs = 0;
  uint64_t total_num_output_subgraphs = 0;
  double total_time_in_milliseconds = 0;

  // to be updated in derived class
  uint64_t intersection_count = 0;
  uint64_t total_intersection_input_size = 0;
  uint64_t total_intersection_output_size = 0;
  uint64_t distinct_intersection_count =
      0;  // the minimal number of intersection needed if all intersection function call results can be cached

  inline CompressedSubgraphs const& getCurrentInput() {return (*current_inputs_)[input_index_];}
  inline CompressedSubgraphs const& getPreviousInput() {return (*current_inputs_)[input_index_-1];}
  inline bool hasNextInput() {return input_index_ < input_end_index_;}
  inline uint32_t& getInputIndex() {input_index_;}
  inline setInputIndex(uint32_t input_index) {input_index_ = input_index;}
  inline setInputEndIndex(uint32_t input_end_index) {input_end_index_ = input_end_index;}
  inline setCurrentInputs(td::vector<CompressedSubgraphs>* current_inputs) {current_inputs_ = current_inputs;}
  inline uint32_t getTotalInputSize() {return total_input_size_ - (current_inputs_size - input_index_);}
  
  private:
    uint32_t input_index_;
    uint32_t input_end_index_;
    const std::vector<CompressedSubgraphs>* current_inputs_ = nullptr;
}

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

  inline virtual void setCandidateSets(const CandidateSetView* candidates) { candidates_ = candidates; }
  inline const CandidateSetView* getCandidateSet() const { return candidates_; }
  inline const uint32_t getInputIndex() const { return input_index_; }
  inline const auto& getSameLabelKeyIndices() const { return same_label_key_indices_; }
  inline const auto& getSameLabelSetIndices() const { return same_label_set_indices_; }
  inline auto getSetPruningThreshold() const { return set_pruning_threshold_; }
  inline QueryVertexID getTargetQueryVertex() const { return target_vertex_; }

  inline const auto& getMatchingOrderIndices() const { return matching_order_indices_; }
  inline void setMatchingOrderIndices(std::vector<std::pair<bool, uint32_t>>&& matching_order_indices) {
    matching_order_indices_ = std::move(matching_order_indices);
  }

  inline bool filter(const CompressedSubgraphs& subgraphs) {
    DCHECK(subgraph_filter_ != nullptr);
    return subgraph_filter_->filter(subgraphs);
  }

  inline uint32_t filter(std::vector<CompressedSubgraphs>& subgraphs, uint32_t start, uint32_t end) {
    DCHECK(subgraph_filter_ != nullptr);
    return subgraph_filter_->filter(subgraphs, start, end);
  }

  virtual std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<CandidateSetView>& candidate_sets) = 0;

  // for backward compabilitity of tests
  inline std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<std::vector<VertexID>>& candidate_sets) {
    std::vector<CandidateSetView> views(candidate_sets.begin(), candidate_sets.end());
    return computeBipartiteGraphs(g, views);
  }

  virtual std::vector<std::unique_ptr<GraphPartitionBase>> computeGraphPartitions(
      const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const = 0;

  virtual void input(const std::vector<CompressedSubgraphs>& inputs, uint32_t input_index, 
                      uint32_t input_end_index, const void* data_graph, TraverseContext* ctx) const {
    ctx->setCurrentInputs(&inputs);
    ctx->setInputIndex(input_index);
    ctx->setInputEndIndex(input_end_index);
    ctx->setCurrentInputs(data_graph); 
  }

  virtual uint32_t expand(uint32_t cap, TraverseContext* ctx) const = 0;

  virtual void inputAndProfile(const std::vector<CompressedSubgraphs>& inputs, uint32_t input_index, 
                      uint32_t input_end_index, const void* data_graph, TraverseContext* ctx) const {
    input(inputs, input_index, input_end_index, data_graph, ctx);
    ctx->current_inputs_size = (input_end_index - input_index));
    ctx->total_input_size += ctx->current_inputs_size;
  }

  uint32_t expandAndProfile(uint32_t cap, uint32_t query_type, TraverseContext* ctx) const {
    auto start = std::chrono::high_resolution_clock::now();
    auto n = expandAndProfileInner(cap, query_type, ctx);
    auto stop = std::chrono::high_resolution_clock::now();
    total_time_in_milliseconds_ +=
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

  void updateIntersectInfo(uint32_t input_size, uint32_t output_size, TraverseContext* ctx) const {
    ++ctx->intersection_count;
    ctx->total_intersection_input_size += input_size;
    ctx->total_intersection_output_size += output_size;
  }

  std::string toProfileString(TraverseContext* ctx) const override {
    std::stringstream ss;
    ss << toString() << ',' << ctx->total_time_in_milliseconds << ',' << getTotalInputSize() << ',' << ctx->total_output_size
       << ',' << ctx->total_num_input_subgraphs << ',' << ctx->total_num_output_subgraphs << ',' << ctx->intersection_count << ','
       << ctx->total_intersection_input_size << ',' << ctx->total_intersection_output_size << ','
       << ctx->distinct_intersection_count;
    return ss.str();
  }

 protected:
  virtual uint32_t expandAndProfileInner(uint32_t cap, uint32_t query_type, TraverseContext* ctx) const = 0;
  inline uint32_t getTotalInputSize(TraverseContext* ctx) const {return ctx->getTotalInputSize();}
};

}  // namespace circinus

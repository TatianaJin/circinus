// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
// the License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/traverse_operator.h"
#include "utils/flags.h"

namespace circinus {

class SplitSetContext : public TraverseContext {
 private:
  uint32_t set_offset_ = 0;

 public:
  SplitSetContext(const void* data_graph, std::vector<CompressedSubgraphs>* outputs, QueryType type)
      : TraverseContext(data_graph, outputs, type) {}

  void setSplitOffset(uint32_t offset) { set_offset_ = offset; }
  uint32_t getSplitOffset() const { return set_offset_; }

  std::unique_ptr<TraverseContext> clone() const override { return std::make_unique<SplitSetContext>(*this); }
};

// FIXME(tatiana): for now use the traverse operator interface while this op does not really traverse, operator should
// have an interface for processing CompressedSubgraphs
class SplitSetOperator : public TraverseOperator {
 private:
  uint32_t split_set_index_;
  uint32_t split_size_;

 public:
  // TODO(tatiana): now directly using flag
  SplitSetOperator(uint32_t split_set_index, uint32_t root_order)
      : TraverseOperator(DUMMY_QUERY_VERTEX, nullptr),
        split_set_index_(split_set_index),
        split_size_(std::pow(FLAGS_batch_size, 1. / root_order)) {
    split_size_ = std::max(split_size_, 1u);
  }

  bool extend_vertex() const override { return false; }

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    if (ctx->getQueryType() == QueryType::Profile) return expandInner<QueryType::Profile>(batch_size, ctx);
    if (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect)
      return expandInner<QueryType::ProfileCandidateSIEffect>(batch_size, ctx);
    CHECK(ctx->getQueryType() == QueryType::ProfileWithMiniIntersection) << "unknown query type "
                                                                         << (uint32_t)ctx->getQueryType();
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, ctx);
  }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    return std::make_unique<SplitSetContext>(graph, outputs, profile);
  }

  std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<CandidateSetView>& candidate_sets) override {
    return std::vector<std::unique_ptr<BipartiteGraph>>();
  }

  std::vector<std::unique_ptr<GraphPartitionBase>> computeGraphPartitions(
      const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const override {
    return std::vector<std::unique_ptr<GraphPartitionBase>>();
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "SplitSetOperator (split set index " << split_set_index_ << ", split size " << split_size_ << ')';
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first, input_key_size.second};
  }

  void setPartialOrder(const PartialOrder& po, const unordered_map<QueryVertexID, uint32_t>& seen_vertices) override {
    // do nothing
  }

 private:
  template <QueryType profile>
  uint32_t expandInner(uint32_t batch_size, TraverseContext* base_ctx) const {
    uint32_t n_outputs = 0;
    auto ctx = dynamic_cast<SplitSetContext*>(base_ctx);
    DCHECK(ctx != nullptr) << "Expect pointer to SplitSetContext but got " << getTypename(*ctx);

    while (ctx->hasNextInput() && n_outputs < batch_size) {
      const auto& input = ctx->getCurrentInput();
      DCHECK_LT(split_set_index_, input.getNumSets());
      auto& set = input.getSet(split_set_index_);

      uint32_t index = ctx->getSplitOffset();
      if (index >= set->size()) {
        ctx->nextInput();
        ctx->setSplitOffset(0);
        continue;
      }
      while (index + split_size_ <= set->size()) {
        auto& output = ctx->copyOutput(input);
        output.UpdateSets(split_set_index_, VertexSet(set, index, split_size_));
        index += split_size_;
        ++n_outputs;
        if (n_outputs == batch_size) {
          ctx->setSplitOffset(index);
          return n_outputs;
        }
      }
      ctx->copyOutput(input).UpdateSets(split_set_index_, VertexSet(set, index, set->size() - index));
      ++n_outputs;
      ctx->nextInput();
      ctx->setSplitOffset(0);
    }
    return n_outputs;
  }
};

}  // namespace circinus

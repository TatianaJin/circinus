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
#include <utility>
#include <vector>

#include "graph/candidate_set_view.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/partitioned_graph.h"
#include "ops/filters/filter.h"
#include "ops/operator.h"

namespace circinus {

class InputOperator : public Operator {
 protected:
  const QueryVertexID starting_vertex_;
  const bool inputs_are_keys_;

 public:
  InputOperator(const QueryVertexID starting_vertex, const bool inputs_are_keys)
      : starting_vertex_(starting_vertex), inputs_are_keys_(inputs_are_keys) {}

  virtual std::vector<CompressedSubgraphs> getInputs(const void* g, const std::vector<CandidateSetView>& candidates) {
    if (inputs_are_keys_) {
      return std::vector<CompressedSubgraphs>(candidates[starting_vertex_].begin(), candidates[starting_vertex_].end());
    }
    return std::vector<CompressedSubgraphs>({CompressedSubgraphs(std::make_shared<std::vector<VertexID>>(
        std::vector<VertexID>(candidates[starting_vertex_].begin(), candidates[starting_vertex_].end())))});
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "InputOperator: query vertex " << starting_vertex_ << " key " << inputs_are_keys_;
    return ss.str();
  }
};

class PartitionedInputOperator : public InputOperator {
 private:
  const std::vector<std::pair<QueryVertexID, QueryVertexID>>* qv_pivots_;

 public:
  PartitionedInputOperator(const QueryVertexID starting_vertex, const bool inputs_are_keys,
                           const std::vector<std::pair<QueryVertexID, QueryVertexID>>* qv_pivots)
      : InputOperator(starting_vertex, inputs_are_keys), qv_pivots_(qv_pivots) {}

  std::vector<CompressedSubgraphs> getInputs(const void* data_graph,
                                             const std::vector<CandidateSetView>& candidates) override {
    if (!qv_pivots_->empty()) {
      auto g = reinterpret_cast<const ReorderedPartitionedGraph*>(data_graph);
      std::vector<VertexID> candidate[2];
      uint32_t cur_idx = 0;
      for (uint32_t i = 1; i < qv_pivots_->size(); ++i) {
        CHECK_EQ((*qv_pivots_)[i].second, (*qv_pivots_)[i - 1].first) << "qv pivots is not a path.";
      }

      candidate[cur_idx].assign(candidates[(*qv_pivots_)[cur_idx].second].begin(),
                                candidates[(*qv_pivots_)[cur_idx].second].end());
      for (auto& qv_pivot_pair : *qv_pivots_) {
        for (auto data_vertex : candidates[qv_pivot_pair.first]) {
          const auto& nbrs = g->getOutNeighbors(data_vertex);
          candidate[cur_idx ^ 1].clear();
          if (NeighborhoodFilter::intersectionNotNull(
                  std::make_pair(candidate[cur_idx].data(), candidate[cur_idx].size()), nbrs)) {
            candidate[cur_idx ^ 1].push_back(data_vertex);
          }
          cur_idx ^= 1;
        }
      }
      if (inputs_are_keys_) {
        return std::vector<CompressedSubgraphs>(candidate[cur_idx].begin(), candidate[cur_idx].end());
      } else {
        return std::vector<CompressedSubgraphs>(
            {CompressedSubgraphs(std::make_shared<std::vector<VertexID>>(std::move(candidate[cur_idx])))});
      }
    }
    return InputOperator::getInputs(data_graph, candidates);
  }
};

}  // namespace circinus

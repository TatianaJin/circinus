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
#include "ops/filters/filter.h"
#include "ops/operator.h"

namespace circinus {

class InputOperator : public Operator {
 private:
  std::vector<std::pair<QueryVertexID, QueryVertexID>>* qv_pivots_;
  QueryVertexID starting_vertex_;
  const bool inputs_are_keys_;

 public:
  InputOperator(const QueryVertexID starting_vertex, const bool inputs_are_keys,
                std::vector<std::pair<QueryVertexID, QueryVertexID>>* qv_pivots)
      : starting_vertex_(starting_vertex), inputs_are_keys_(inputs_are_keys), qv_pivots_(qv_pivots) {}

  virtual std::vector<CompressedSubgraphs> getInputs(const Graph* data_graph,
                                                     const std::vector<CandidateSetView>& candidates) {
    std::vector<VertexID> candidate[2];
    uint32_t cur_idx = 0;
    if (qv_pivots_->size() != 0) {
      candidate[0].reserve(candidates[(*qv_pivots_)[0].second].size());
      for (auto query_vertex : candidates[(*qv_pivots_)[0].second]) {
        candidate[0].emplace_back(query_vertex);
      }
      for (auto& qv_pivot_pair : *qv_pivots_) {
        for (auto query_vertex : candidates[qv_pivot_pair.first]) {
          const auto& nbrs = data_graph->getOutNeighbors(query_vertex);
          candidate[cur_idx ^ 1].clear();
          if (NeighborhoodFilter::intersectionNotNull(
                  std::make_pair(candidate[cur_idx].data(), candidate[cur_idx].size()), nbrs)) {
            candidate[cur_idx ^ 1].push_back(query_vertex);
          }
          cur_idx ^= 1;
        }
      }
    } else {
      for (auto query_vertex : candidates[starting_vertex_]) {
        candidate[cur_idx].push_back(query_vertex);
      }
    }

    if (inputs_are_keys_) {
      return std::vector<CompressedSubgraphs>(candidate[cur_idx].begin(), candidate[cur_idx].end());
    } else {
      return std::vector<CompressedSubgraphs>(
          {CompressedSubgraphs(std::make_shared<std::vector<VertexID>>(candidate[cur_idx]))});
    }
  }
};

}  // namespace circinus

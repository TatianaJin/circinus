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

  virtual std::vector<CompressedSubgraphs> getInputs(const void* g,
                                                     const std::vector<CandidateSetView>& candidates) const;

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
                                             const std::vector<CandidateSetView>& candidates) const override;

  std::string toString() const override {
    std::stringstream ss;
    ss << "PartitionedInputOperator: query vertex " << starting_vertex_ << " key " << inputs_are_keys_ << " n_pivots "
       << qv_pivots_->size();
    return ss.str();
  }
};

}  // namespace circinus

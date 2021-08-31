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
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "graph/candidate_set_view.h"
#include "graph/compressed_subgraphs.h"
#include "ops/operator.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

class InputOperator : public Operator {
 protected:
  const QueryVertexID starting_vertex_;
  const bool inputs_are_keys_;

 public:
  InputOperator(const QueryVertexID starting_vertex, const bool inputs_are_keys)
      : starting_vertex_(starting_vertex), inputs_are_keys_(inputs_are_keys) {}

  virtual ~InputOperator() {}

  inline std::pair<uint32_t, uint32_t> getOutputSize() const { return {inputs_are_keys_, 1}; }

  virtual std::vector<CompressedSubgraphs> getInputs(const void* g,
                                                     const std::vector<CandidateSetView>& candidates) const;

  // TODO(profile): record intersection count in input op?
  virtual void inputAndProfile(const void* g, const std::vector<CandidateSetView>& candidates,
                               std::vector<CompressedSubgraphs>* output, ProfileInfo* info) const {
    auto start = std::chrono::high_resolution_clock::now();
    *output = getInputs(g, candidates);
    LOG(INFO) << "------------";
    auto end = std::chrono::high_resolution_clock::now();
    info->total_time_in_milliseconds += toMilliseconds(start, end);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "InputOperator: query vertex " << starting_vertex_ << " key " << inputs_are_keys_;
    return ss.str();
  }

  std::string toProfileString(const ProfileInfo& info) const override {
    std::stringstream ss;
    ss << toString() << ',' << info.total_time_in_milliseconds;
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

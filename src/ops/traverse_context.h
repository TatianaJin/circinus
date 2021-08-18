// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License");  you may not use this file except in compliance with
// the License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions andlimitations under the License.

#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "graph/candidate_set_view.h"
#include "graph/compressed_subgraphs.h"
#include "graph/types.h"
#include "utils/query_utils.h"

namespace circinus {

class TraverseContext : public ProfileInfo {
 private:
  uint32_t input_index_ = 0;
  uint32_t input_start_index_ = 0;
  uint32_t input_end_index_ = 0;
  const std::vector<CompressedSubgraphs>* current_inputs_ = nullptr;

  uint32_t output_idx_ = 0;
  std::vector<CompressedSubgraphs>* outputs_ = nullptr;

  const void* current_data_graph_ = nullptr;
  QueryType query_type_ = QueryType::Execute;

 public:
  TraverseContext() {}
  TraverseContext(const void* data_graph, std::vector<CompressedSubgraphs>* outputs, QueryType type)
      : outputs_(outputs), current_data_graph_(data_graph), query_type_(type) {}

  virtual ~TraverseContext() {}

  virtual std::unique_ptr<TraverseContext> clone() const = 0;

  inline const CompressedSubgraphs& getCurrentInput() const { return (*current_inputs_)[input_index_]; }
  inline const CompressedSubgraphs& getPreviousInput() const { return (*current_inputs_)[input_index_ - 1]; }
  inline bool hasNextInput() const { return input_index_ < input_end_index_; }
  inline uint32_t getInputIndex() const { return input_index_; }
  inline uint32_t getTotalInputSize() const { return input_index_ - input_start_index_; }
  inline auto getOutputs() const { return outputs_; }
  inline auto getOutputSize() const { return output_idx_; }
  inline auto getQueryType() const { return query_type_; }

  template <typename G>
  inline const G* getDataGraph() const {
    return reinterpret_cast<const G*>(current_data_graph_);
  }

  inline void nextInput() { ++input_index_; }

  inline void setInput(const std::vector<CompressedSubgraphs>& inputs, uint32_t start, uint32_t end) {
    current_inputs_ = &inputs;
    input_start_index_ = input_index_ = start;
    input_end_index_ = end;
  }

  inline void setOutputBuffer(std::vector<CompressedSubgraphs>& outputs) { outputs_ = &outputs; }
  inline void popOutputs(uint32_t size) { output_idx_ -= size; }
  inline void popOutput() { --output_idx_; }
  inline void resetOutputs() { output_idx_ = 0; }

  template <typename... Args>
  inline CompressedSubgraphs* newOutput(Args&&... args) {
    DCHECK_LT(output_idx_, outputs_->size()) << "output buffer addr " << outputs_;
    auto ret = (*outputs_)[output_idx_].reset(std::forward<Args>(args)...);
    // do not increment output_idx_ if the current output is invalid
    output_idx_ += (ret != nullptr);
    return ret;
  }

  inline CompressedSubgraphs& copyOutput(const CompressedSubgraphs& from) {
    (*outputs_)[output_idx_] = from;
    return (*outputs_)[output_idx_++];
  }
};

}  // namespace circinus

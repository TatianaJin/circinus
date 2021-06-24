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

#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/graph_view.h"
#include "ops/operator.h"
#include "ops/operators.h"
#include "utils/flags.h"
#include "utils/profiler.h"

namespace circinus {

// forward declaration
class Task;

/** Uses flag: FLAGS_batch_size */
class OperatorTree {
  std::vector<Operator*> operators_;  // FIXME(tatiana): use unique_ptr?
  Profiler* profiler_;

 public:
  ~OperatorTree() { clear(); }

  void clear() {
    for (auto op : operators_) {
      delete op;
    }
    operators_.clear();
  }
  inline bool empty() const { return operators_.empty(); }

  inline Operator* root() const { return operators_.front(); }
  inline void push_back(Operator* t) { operators_.push_back(t); }
  inline void reserve(uint32_t size) { operators_.reserve(size); }
  inline void setProfiler(Profiler* profiler) { profiler_ = profiler; }
  inline Operator* getOperator(uint32_t idx) const { return operators_[idx]; }
  inline size_t getOperatorSize() const { return operators_.size(); }

  inline void setOutput(Outputs* outputs) const {
    auto op = operators_.back();
    auto output_op = dynamic_cast<OutputOperator*>(op);
    CHECK(output_op != nullptr) << op->toString();
    output_op->setOutput(outputs);
  }

  bool execute(const Graph* g, const std::vector<CompressedSubgraphs>& inputs, uint32_t level = 0);
  bool profile(const Graph* g, const std::vector<CompressedSubgraphs>& inputs, uint32_t query_type, uint32_t level = 0);

  template <typename G>
  bool execute(const std::vector<GraphView<G>*>& data_graphs_for_operators,
               const std::vector<CompressedSubgraphs>& inputs, uint32_t level = 0) {
    std::vector<CompressedSubgraphs> outputs;
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      return output_op->validateAndOutput(inputs, 0);
    }
    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    traverse_op->input(inputs, data_graphs_for_operators[level]);
    while (true) {
      outputs.clear();
      auto size = traverse_op->expand(&outputs, FLAGS_batch_size);
      if (size == 0) {
        break;
      }
      if (execute<G>(data_graphs_for_operators, outputs, level + 1)) {
        return true;
      }
    }
    return false;
  }

  template <typename G>
  bool profile(const std::vector<GraphView<G>*>& data_graphs_for_operators,
               const std::vector<CompressedSubgraphs>& inputs, uint32_t query_type, uint32_t level = 0) {
    std::vector<CompressedSubgraphs> outputs;
    auto op = operators_[level];
    if (level == operators_.size() - 1) {
      auto output_op = dynamic_cast<OutputOperator*>(op);
      return output_op->validateAndOutputAndProfile(inputs, 0);
    }
    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    traverse_op->inputAndProfile(inputs, data_graphs_for_operators[level]);
    while (true) {
      outputs.clear();
      auto size = traverse_op->expandAndProfile(&outputs, FLAGS_batch_size, query_type);
      if (size == 0) {
        break;
      }
      if (profile<G>(data_graphs_for_operators, outputs, query_type, level + 1)) {
        return true;
      }
    }
    return false;
  }
};

}  // namespace circinus

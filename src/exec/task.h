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

namespace circinus {

// TODO(tatiana): add query id for concurrent processing
class Task {
  uint32_t level_;
  std::vector<CompressedSubgraphs> input_;
  const Graph* data_graph_;

 public:
  Task(uint32_t level, std::vector<CompressedSubgraphs>&& input, const Graph* graph)
      : level_(level), input_(std::move(input)), data_graph_(graph) {}

  inline uint32_t getLevel() const { return level_; }
  inline const auto& getInput() const { return input_; }
  inline const Graph* getDataGraph() const { return data_graph_; }
};

}  // namespace circinus

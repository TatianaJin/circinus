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

#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/types.h"
#include "ops/operator.h"

namespace circinus {

/** set1 and set2 must be sorted in ascending order */
void intersect(const std::pair<const VertexID*, uint32_t>& set1, const std::pair<const VertexID*, uint32_t>& set2,
               std::vector<VertexID>* intersection);

/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection) {
  intersect(std::make_pair(set1.data(), (uint32_t)set1.size()), set2, intersection);
}

/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::vector<VertexID>& set1, const std::vector<VertexID>& set2,
                      std::vector<VertexID>* intersection) {
  intersect(set1, std::make_pair(set2.data(), (uint32_t)set2.size()), intersection);
}

void intersectInplace(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection);

class TraverseOperator : public Operator {
 protected:
  const std::vector<VertexID>* candidates_;

  /* transient variables for recording the current inputs */
  uint32_t input_index_ = 0;
  const std::vector<CompressedSubgraphs>* current_inputs_;
  const Graph* current_data_graph_;

 public:
  virtual ~TraverseOperator() {}

  inline void setCandidateSets(const std::vector<VertexID>* candidates) { candidates_ = candidates; }

  virtual void input(const std::vector<CompressedSubgraphs>& inputs, const Graph* data_graph) {
    current_inputs_ = &inputs;
    input_index_ = 0;
    current_data_graph_ = data_graph;
  }

  virtual uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) = 0;
};

}  // namespace circinus

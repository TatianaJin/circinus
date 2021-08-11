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
#include <unordered_set>
#include <vector>

#include "ops/traverse_operator.h"
#include "ops/traverse_operator_utils.h"

namespace circinus {

#ifdef INTERSECTION_CACHE
class ExpandVertexTraverseContext : public TraverseContext, public MultiparentIntersectionCache {
#else
class ExpandVertexTraverseContext : public TraverseContext {
#endif

 protected:
  const CandidateSetView* candidates_ = nullptr;

  // for profiling
  std::vector<unordered_set<std::string>> parent_tuple_sets_;

 public:
  ExpandVertexTraverseContext(const CandidateSetView* candidates, const void* graph,
                              std::vector<CompressedSubgraphs>* outputs, QueryType profile, uint32_t parent_size)
      : TraverseContext(graph, outputs, profile), candidates_(candidates), parent_tuple_sets_(parent_size) {}

  virtual ~ExpandVertexTraverseContext() {}

  std::unique_ptr<TraverseContext> clone() const override {
    return std::make_unique<ExpandVertexTraverseContext>(*this);
  }

  inline const CandidateSetView* getCandidateSet() const { return candidates_; }

  inline void updateDistinctSICount(uint32_t depth, std::vector<VertexID>& parent_tuple, uint32_t pidx) {
    distinct_intersection_count +=
        parent_tuple_sets_[depth].emplace((const char*)parent_tuple.data(), (pidx + 1) * sizeof(VertexID)).second;
  }
};

}  // namespace circinus

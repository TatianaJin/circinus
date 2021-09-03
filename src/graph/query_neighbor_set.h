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

#include "graph/query_graph.h"
#include "graph/types.h"

namespace circinus {

class QueryNeighborSet {
  const QueryVertexID* start_ = nullptr;
  const QueryVertexID* end_ = nullptr;

 public:
  using ConstIterator = const QueryVertexID*;
  QueryNeighborSet() {}

  inline void addRange(const QueryVertexID* start, const QueryVertexID* end) {
    DCHECK(start != nullptr);
    DCHECK(end != nullptr);
    CHECK_LT(start, end);
    start_ = start;
    end_ = end;
  }
  inline ConstIterator begin() const { return start_; }
  inline ConstIterator end() const { return end_; }
  inline size_t size() const { return start_ == nullptr ? 0 : std::distance(start_, end_); }
  inline bool empty() const { return start_ == nullptr; }

  friend bool operator<(const QueryNeighborSet& a, const QueryNeighborSet& b) { return a.size() < b.size(); }
};

}  // namespace circinus

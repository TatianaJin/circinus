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

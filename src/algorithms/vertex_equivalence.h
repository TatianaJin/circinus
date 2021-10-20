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

#include <functional>
#include <unordered_set>
#include <utility>

#include "algorithms/partial_order.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

class VertexEquivalence {
#ifdef USE_STL
  struct HashPair {
    template <class T1, class T2>
    inline std::size_t operator()(const std::pair<T1, T2>& pair) const {
      return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
  };

  std::unordered_set<std::pair<QueryVertexID, QueryVertexID>, HashPair> equivalent_pairs_;
#else
  phmap::flat_hash_set<std::pair<QueryVertexID, QueryVertexID> > equivalent_pairs_;
#endif

  const QueryGraph* q_ = nullptr;

 public:
  VertexEquivalence(const QueryGraph& q, PartialOrder* po);

  inline const QueryGraph* getQueryGraph() const { return q_; }

  inline bool empty() const { return equivalent_pairs_.empty(); }

  inline bool isEquivalent(QueryVertexID u1, QueryVertexID u2) const {
    auto pair = (u1 < u2) ? std::make_pair(u1, u2) : std::make_pair(u2, u1);
    return equivalent_pairs_.count(pair) == 1;
  }

 private:
  static bool setsEqualExcept(const unordered_set<QueryVertexID>& set1, const unordered_set<QueryVertexID>& set2,
                              QueryVertexID u1, QueryVertexID u2);
};

}  // namespace circinus

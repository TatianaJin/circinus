#pragma once

#include <algorithm>
#include <unordered_set>
#include <vector>

#include "algorithms/search.h"
#include "graph/candidate_set_view.h"
#include "graph/types.h"
#include "graph/vertex_set_view.h"
#include "utils/hashmap.h"

namespace circinus {

void _leapfrogJoinPointers(std::array<std::pair<const VertexID*, const VertexID*>, 2>& pointers,
                           std::vector<VertexID>* intersection, const unordered_set<VertexID>& except);

inline void leapfrogJoin(const SingleRangeVertexSetView& set1, const CandidateSetView& set2,
                         std::vector<VertexID>* intersection, const unordered_set<VertexID>& except) {
  if (set1.empty() || set2.empty() || set1.front() > set2.back() || set2.front() > set1.back()) return;

  std::array<std::pair<const VertexID*, const VertexID*>, 2> pointers;
  pointers[0] = {set1.begin(), set1.end()};
  pointers[1] = {set2.begin(), set2.end()};
  _leapfrogJoinPointers(pointers, intersection, except);
}

inline void leapfrogJoin(const SingleRangeVertexSetView& set1, const SingleRangeVertexSetView& set2,
                         std::vector<VertexID>* intersection, const unordered_set<VertexID>& except) {
  if (set1.empty() || set2.empty() || set1.front() > set2.back() || set2.front() > set1.back()) return;

  std::array<std::pair<const VertexID*, const VertexID*>, 2> pointers;
  pointers[0] = {set1.begin(), set1.end()};
  pointers[1] = {set2.begin(), set2.end()};
  _leapfrogJoinPointers(pointers, intersection, except);
}

void leapfrogJoin(std::vector<SingleRangeVertexSetView>& sets, std::vector<VertexID>* intersection,
                  const unordered_set<VertexID>& except);

void leapfrogJoin(std::vector<VertexSetView>& sets, std::vector<VertexID>* intersection,
                  const unordered_set<VertexID>& except);

void leapfrogJoin(const SingleRangeVertexSetView& set1, const VertexSetView& set2, std::vector<VertexID>* intersection,
                  const unordered_set<VertexID>& except);

}  // namespace circinus

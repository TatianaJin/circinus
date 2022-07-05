#include "algorithms/leapfrog_join.h"

#include <algorithm>
#include <unordered_set>
#include <vector>

namespace circinus {

void _leapfrogJoinPointers(std::array<std::pair<const VertexID*, const VertexID*>, 2>& pointers,
                           std::vector<VertexID>* intersection, const unordered_set<VertexID>& except) {
  uint32_t p = *pointers[0].first > *pointers[1].first;  // index to smaller pointer
  VertexID x_prime, x;
  while (true) {
    x_prime = *pointers[1 - p].first;
    x = *pointers[p].first;
    if (x == x_prime) {  // match
      if (except.count(x) == 0) {
        intersection->push_back(x);
      }
      ++pointers[p].first;  // next
    } else {                // x < x_prime, seek set p to find x_prime or a bit larger
      pointers[p].first = lowerBound(pointers[p].first, pointers[p].second, x_prime);
    }
    if (pointers[p].first == pointers[p].second) {
      break;  // at end
    }
    p = 1 - p;  // now the smaller pointer changes
  }
}

void leapfrogJoin(std::vector<SingleRangeVertexSetView>& sets, std::vector<VertexID>* intersection,
                  const unordered_set<VertexID>& except) {
  for (auto& set : sets) {
    if (set.empty()) {
      return;
    }
  }

  std::sort(sets.begin(), sets.end(),
            [](const SingleRangeVertexSetView& a, const SingleRangeVertexSetView& b) { return a.front() < b.front(); });
  std::vector<typename SingleRangeVertexSetView::ConstIterator> iters;
  iters.reserve(sets.size());
  for (auto& set : sets) {
    iters.push_back(set.begin());
  }

  uint32_t p = 0;
  uint32_t n_sets = iters.size();
  VertexID x_prime, x;
  while (true) {
    x_prime = *iters[(p + n_sets - 1) % n_sets];
    x = *iters[p];
    if (x == x_prime) {  // match
      if (except.count(x) == 0) {
        intersection->push_back(x);
      }
      ++iters[p];  // next
    } else {       // x < x_prime, seek iter p
      iters[p] = lowerBound(iters[p], sets[p].end(), x_prime);
    }
    if (iters[p] == sets[p].end()) {
      break;
    }
    p = (p + 1) % n_sets;
  }
}

void leapfrogJoin(std::vector<VertexSetView>& sets, std::vector<VertexID>* intersection,
                  const unordered_set<VertexID>& except) {
  bool all_single_range = true;
  for (auto& set : sets) {
    if (set.empty()) return;
    if (set.getRanges().size() != 1) {
      all_single_range = false;
    }
  }
  if (all_single_range) {
    std::vector<SingleRangeVertexSetView> ssets;
    ssets.reserve(sets.size());
    for (auto& set : sets) {
      ssets.emplace_back(set.getRanges().front().first, set.getRanges().front().second);
    }
    return leapfrogJoin(ssets, intersection, except);
  }

  std::sort(sets.begin(), sets.end(),
            [](const VertexSetView& a, const VertexSetView& b) { return a.front() < b.front(); });
  std::vector<typename VertexSetView::ConstIterator> iters;
  iters.reserve(sets.size());
  for (auto& set : sets) {
    iters.push_back(set.begin());
  }

  uint32_t p = 0;
  uint32_t n_sets = iters.size();
  VertexID x_prime, x;
  while (true) {
    x_prime = *iters[(p + n_sets - 1) % n_sets];
    x = *iters[p];
    if (x == x_prime) {  // match
      if (except.count(x) == 0) {
        intersection->push_back(x);
      }
      ++iters[p];  // next
    } else {       // x < x_prime, seek iter p
      iters[p].seek(x_prime);
    }
    if (iters[p] == sets[p].end()) {
      break;
    }
    p = (p + 1) % n_sets;
  }
}

void leapfrogJoin(const SingleRangeVertexSetView& set1, const VertexSetView& set2, std::vector<VertexID>* intersection,
                  const unordered_set<VertexID>& except) {
  if (set2.getRanges().size() == 1) {
    return leapfrogJoin(set1, SingleRangeVertexSetView(set2.getRanges().front().first, set2.getRanges().front().second),
                        intersection, except);
  }
  std::vector<VertexSetView> sets{set1, set2};
  leapfrogJoin(sets, intersection, except);
}

}  // namespace circinus

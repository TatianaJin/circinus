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

#include "ops/filters/subgraph_filter.h"

#include <algorithm>
#include <vector>

#include "graph/compressed_subgraphs.h"

namespace circinus {

class SetPrunningSubgraphFilter : public SubgraphFilter {
  const std::vector<std::vector<uint32_t>> pruning_set_indices_;
  std::vector<uint64_t> set_pruning_thresholds_;

 public:
  explicit SetPrunningSubgraphFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices)
      : pruning_set_indices_(std::move(pruning_set_indices)) {}

  void setPruningSetThresholds(std::vector<uint64_t>&& pruning_set_thresholds) override {
    set_pruning_thresholds_ = std::move(pruning_set_thresholds);
    CHECK_EQ(set_pruning_thresholds_.size(), pruning_set_indices_.size());
  }

  const std::vector<uint32_t>* getPruningSets(uint32_t idx) const override {
    DCHECK_LT(idx, pruning_set_indices_.size());
    return &pruning_set_indices_[idx];
  }

  uint64_t getSetPruningThreshold(uint32_t idx) const override {
    DCHECK_LT(idx, set_pruning_thresholds_.size());
    return set_pruning_thresholds_[idx];
  }

  uint32_t filter(std::vector<CompressedSubgraphs>& subgraphs, uint32_t start, uint32_t end) override {
    auto end_iter = subgraphs.begin() + end;
    auto start_iter = subgraphs.begin() + start;
    auto old_size = subgraphs.size();
    subgraphs.erase(std::remove_if(start_iter, end_iter, [this](auto& group) { return filter(group); }), end_iter);
    return old_size - subgraphs.size();
  }

  /** @returns True if pruned. */
  bool filter(const CompressedSubgraphs& subgraphs) override {
    for (auto& pruning_sets : pruning_set_indices_) {
      // get pruning set pointers sorted by size in ascending order
      std::vector<std::vector<VertexID>*> pruning_sets_sorted_by_size;
      uint32_t n_sets = pruning_sets.size();
      pruning_sets_sorted_by_size.reserve(n_sets);
      for (auto index : pruning_sets) {
        auto ptr = subgraphs.getSet(index).get();
        // only check sets whose size is smaller than the number of pruning sets
        if (ptr->size() < n_sets) pruning_sets_sorted_by_size.push_back(ptr);
      }
      if (pruning_sets_sorted_by_size.empty()) {
        continue;
      }
      std::sort(pruning_sets_sorted_by_size.begin(), pruning_sets_sorted_by_size.end(),
                [](const auto& set1, const auto& set2) { return set1->size() < set2->size(); });
      // dfs to check whether there is at least one tuple without repeated vertex
      if (notOneValidTuple(subgraphs, pruning_sets_sorted_by_size)) {
        return true;
      }
    }
    return false;
  }

 private:
  static inline bool notOneValidTuple(const CompressedSubgraphs& subgraphs, std::vector<std::vector<VertexID>*> sets) {
    uint32_t n_sets = sets.size();
    uint32_t last_depth = n_sets - 1;
    uint32_t current_depth = 0;
    std::vector<uint32_t> set_index(n_sets, 0);
    unordered_set<VertexID> existing_vertices;
    existing_vertices.reserve(n_sets - 1);
    while (true) {
      while (set_index[current_depth] < sets[current_depth]->size()) {
        auto v = (*sets[current_depth])[set_index[current_depth]];
        if (existing_vertices.count(v) == 0) {  // v is valid
          if (current_depth == last_depth) {    // found a valid `n_sets`-tuple
            return false;
          }
          existing_vertices.insert(v);
          ++current_depth;
          set_index[current_depth] = 0;  // start from the first vertex in the next set
        } else {
          ++set_index[current_depth];
        }
      }
      if (current_depth == 0) {
        break;
      }
      --current_depth;
      existing_vertices.erase((*sets[current_depth])[set_index[current_depth]]);
      ++set_index[current_depth];
    }
    return true;
  }
};

SubgraphFilter* SubgraphFilter::newSetPrunningSubgraphFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices) {
  return new SetPrunningSubgraphFilter(std::move(pruning_set_indices));
}

SubgraphFilter* SubgraphFilter::newSetPrunningSubgraphFilter(const std::vector<uint32_t>& pruning_set_indices) {
  if (pruning_set_indices.size() < 2) {
    return new SubgraphFilter();
  }
  std::vector<std::vector<uint32_t>> pruning_sets(1);
  pruning_sets.back() = pruning_set_indices;
  return new SetPrunningSubgraphFilter(std::move(pruning_sets));
}

}  // namespace circinus

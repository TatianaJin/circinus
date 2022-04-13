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
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "algorithms/search.h"
#include "graph/compressed_subgraphs.h"
#include "graph/types.h"
#include "utils/utils.h"

namespace circinus {

class SetPrunningSubgraphFilter : public SubgraphFilter {
 protected:
  std::vector<std::vector<uint32_t>> pruning_set_indices_;
  std::vector<uint64_t> set_pruning_thresholds_;

 public:
  SetPrunningSubgraphFilter() {}

  explicit SetPrunningSubgraphFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices)
      : pruning_set_indices_(std::move(pruning_set_indices)) {}

  SetPrunningSubgraphFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices,
                            std::vector<uint64_t>&& pruning_set_thresholds)
      : pruning_set_indices_(std::move(pruning_set_indices)),
        set_pruning_thresholds_(std::move(pruning_set_thresholds)) {
    if (!set_pruning_thresholds_.empty()) CHECK_EQ(set_pruning_thresholds_.size(), pruning_set_indices_.size());
  }

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

  std::unique_ptr<SubgraphFilter> addPartialOrderSubgraphFilter(std::vector<uint32_t>&& lt_conditions,
                                                                std::vector<uint32_t>&& gt_conditions,
                                                                const std::pair<bool, uint32_t>& target) override;

  uint32_t filter(std::vector<CompressedSubgraphs>& subgraphs, uint32_t start, uint32_t end) override {
    auto end_iter = subgraphs.begin() + end;
    auto start_iter = subgraphs.begin() + start;
    start_iter = std::find_if(start_iter, end_iter, [this](auto& group) { return filter(group); });
    if (start_iter != end_iter) {
      for (auto i = start_iter; ++i != end_iter;) {
        if (!filter(*i)) {
          start_iter->swap(*i);
          ++start_iter;
        }
      }
    }
    return std::distance(start_iter, end_iter);
  }

  /** @returns True if pruned. */
  bool filter(CompressedSubgraphs& subgraphs) const override {
    for (auto& pruning_sets : pruning_set_indices_) {
      // get pruning set pointers sorted by size in ascending order
      std::vector<std::vector<VertexID>*> pruning_sets_sorted_by_size;
      uint32_t n_sets = pruning_sets.size();
      if (n_sets == 1) continue;
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
};  // class SetPrunningSubgraphFilter

/** Prune sets according to less/greater than conditions with the newly matched target */
class PartialOrderSubgraphFilter : public SetPrunningSubgraphFilter {
  // set indices
  std::vector<uint32_t> less_than_indices_;
  std::vector<uint32_t> greater_than_indices_;
  std::pair<bool, uint32_t> target_index_;

 public:
  /**
   * @param conditions A vector of inequality conditions, represented as {set index, is_less_than} pairs.
   * @param target A pair {is_key, index}.
   */
  PartialOrderSubgraphFilter(std::vector<uint32_t>&& lt_conditions, std::vector<uint32_t>&& gt_conditions,
                             const std::pair<bool, uint32_t>& target)
      : less_than_indices_(std::move(lt_conditions)),
        greater_than_indices_(std::move(gt_conditions)),
        target_index_(target) {
    LOG(INFO) << "target index = " << target_index_.second << ", is key = " << target_index_.first << " less than "
              << toString(less_than_indices_) << " greater than " << toString(greater_than_indices_);
  }

  PartialOrderSubgraphFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices,
                             std::vector<uint64_t>&& pruning_set_thresholds, std::vector<uint32_t>&& lt_conditions,
                             std::vector<uint32_t>&& gt_conditions, const std::pair<bool, uint32_t>& target)
      : SetPrunningSubgraphFilter(std::move(pruning_set_indices), std::move(pruning_set_thresholds)),
        less_than_indices_(std::move(lt_conditions)),
        greater_than_indices_(std::move(gt_conditions)),
        target_index_(target) {
    LOG(INFO) << "target index = " << target_index_.second << ", is key = " << target_index_.first << " less than "
              << toString(less_than_indices_) << " greater than " << toString(greater_than_indices_);
  }

  /** @returns True if pruned. */
  bool filter(CompressedSubgraphs& subgraphs) const override {
    // apply set pruning first to avoid copy due to partial pruning by partial order
    if (SetPrunningSubgraphFilter::filter(subgraphs)) return true;
    std::vector<uint32_t> pruned_set_indices;
    if (target_index_.first) {  // target is key
      if (pruneByTarget(subgraphs, subgraphs.getKeyVal(target_index_.second), &pruned_set_indices)) {
        return true;
      }
    } else {  // target is set, prune if set has only one element
      auto target_set = subgraphs.getSet(target_index_.second).get();
      if (target_set->size() == 1 && pruneByTarget(subgraphs, target_set->front(), &pruned_set_indices)) {
        return true;
      }
    }
    if (pruned_set_indices.empty()) return false;
    return SetPrunningSubgraphFilter::filter(subgraphs);
  }

 private:
  /**
   * @param pruned_set_indices Outputs which are the indices of sets whose elements are pruned (partially).
   * @returns True if input should be pruned.
   */
  bool pruneByTarget(CompressedSubgraphs& subgraphs, VertexID target, std::vector<uint32_t>* pruned_set_indices) const {
    for (auto index : less_than_indices_) {  // target should be less (smaller) than the valid set elements
      const auto& set = subgraphs.getSet(index);
      auto lb = circinus::lowerBound(set->begin(), set->end(), target);
      if (lb == set->end()) {  // all elements are no larger than target
        return true;
      }
      if (target == *lb) {
        ++lb;
      } else if (lb == set->begin()) {  // no element is to be pruned
        continue;
      }
      // create a new set containing only lb and afterwards, zero copy
      subgraphs.UpdateSets(index, VertexSet(set, lb - set->begin(), set->end() - lb));
      pruned_set_indices->push_back(index);
    }
    for (auto index : greater_than_indices_) {  // target should be greater than the valid set elements
      const auto& set = subgraphs.getSet(index);
      auto lb = circinus::lowerBound(set->begin(), set->end(), target);
      if (lb == set->begin()) {  // target is no greater than any of the elements
        return true;
      }
      if (lb == set->end()) {  // no element is to be pruned
        continue;
      }
      // create a new set containing only elements before lb, zero copy
      subgraphs.UpdateSets(index, VertexSet(set, 0, lb - set->begin()));
      pruned_set_indices->push_back(index);
    }
    return false;
  }
};  // class PartialOrderSubgraphFilter

std::unique_ptr<SubgraphFilter> SetPrunningSubgraphFilter::addPartialOrderSubgraphFilter(
    std::vector<uint32_t>&& lt_conditions, std::vector<uint32_t>&& gt_conditions,
    const std::pair<bool, uint32_t>& target) {
  return std::make_unique<PartialOrderSubgraphFilter>(std::move(pruning_set_indices_),
                                                      std::move(set_pruning_thresholds_), std::move(lt_conditions),
                                                      std::move(gt_conditions), target);
}

std::unique_ptr<SubgraphFilter> SubgraphFilter::newSetPrunningSubgraphFilter(
    std::vector<std::vector<uint32_t>>&& pruning_set_indices) {
  return std::make_unique<SetPrunningSubgraphFilter>(std::move(pruning_set_indices));
}

std::unique_ptr<SubgraphFilter> SubgraphFilter::newSetPrunningSubgraphFilter(
    const std::vector<uint32_t>& pruning_set_indices) {
  if (pruning_set_indices.size() < 2) {
    return std::make_unique<SubgraphFilter>();
  }
  std::vector<std::vector<uint32_t>> pruning_sets(1);
  pruning_sets.back() = pruning_set_indices;
  return std::make_unique<SetPrunningSubgraphFilter>(std::move(pruning_sets));
}

std::unique_ptr<SubgraphFilter> SubgraphFilter::newPartialOrderSubgraphFilter(std::vector<uint32_t>&& lt_conditions,
                                                                              std::vector<uint32_t>&& gt_conditions,
                                                                              const std::pair<bool, uint32_t>& target) {
  return std::make_unique<PartialOrderSubgraphFilter>(std::move(lt_conditions), std::move(gt_conditions), target);
}

}  // namespace circinus

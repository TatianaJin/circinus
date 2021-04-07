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

namespace circinus {

/**
 * Thread-safe functor-style filters
 */
class SubgraphFilter {
 public:
  /** Returns a pointer to SubgraphFilter, which requires manual delete. The instance does no filtering. */
  static SubgraphFilter* newDummyFilter() { return new SubgraphFilter(); }

  /** Returns a pointer to SubgraphFilter, which requires manual delete. The instance filters {@link
   * CompressedSubgraphs} by each group of pruning sets when n sets have a union of less than n elements */
  static SubgraphFilter* newSetPrunningSubgraphFilter(std::vector<std::vector<uint32_t>>&& pruning_set_indices);
  static SubgraphFilter* newSetPrunningSubgraphFilter(const std::vector<uint32_t>& pruning_set_indices);

  virtual ~SubgraphFilter() {}

  virtual void setPruningSetThresholds(std::vector<uint64_t>&& pruning_set_thresholds) {}
  virtual uint64_t getSetPruningThreshold(uint32_t idx) const { return 0; }

  virtual const std::vector<uint32_t>* getPruningSets(uint32_t idx) const { return nullptr; }

  /**
   * @param subgraphs The compressed groups of subgraphs to filter.
   * @returns The number of CompressedSubgraphs pruned.
   */
  inline uint32_t filter(std::vector<CompressedSubgraphs>& subgraphs) { return filter(subgraphs, 0, subgraphs.size()); }

  /**
   * @param subgraphs The compressed groups of subgraphs to filter with range [start,end).
   * @param start The starting index of compressed groups (inclusive).
   * @param end The ending index of compressed groups (exclusive).
   * @returns The number of CompressedSubgraphs pruned.
   */
  virtual uint32_t filter(std::vector<CompressedSubgraphs>& subgraphs, uint32_t start, uint32_t end) { return 0; }

  /** @returns True if pruned. */
  virtual bool filter(const CompressedSubgraphs& subgraphs) { return false; /* dummy filter: do nothing */ }
};

}  // namespace circinus

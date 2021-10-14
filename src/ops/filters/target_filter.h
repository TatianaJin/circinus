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

#include <memory>
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/types.h"

namespace circinus {

/** Thread-safe functor-style filter */
class TargetFilter {
  std::vector<uint32_t> less_than_indices_;
  std::vector<uint32_t> greater_than_indices_;

 public:
  static std::unique_ptr<TargetFilter> newTargetFilter(std::vector<uint32_t>&& lt_conditions,
                                                       std::vector<uint32_t>&& gt_conditions) {
    return std::make_unique<TargetFilter>(std::move(lt_conditions), std::move(gt_conditions));
  }

  TargetFilter(std::vector<uint32_t>&& lt_conditions, std::vector<uint32_t>&& gt_conditions)
      : less_than_indices_(std::move(lt_conditions)), greater_than_indices_(std::move(gt_conditions)) {}

  void filter(std::vector<VertexID>* targets, const CompressedSubgraphs& group) const {
    // get lt and gt constraints for targets
    VertexID less_than = std::numeric_limits<VertexID>::max();  // min of all constraints
    for (auto index : less_than_indices_) {
      auto val = group.getKeyVal(index);
      less_than = std::min(less_than, val);
    }
    auto size = targets->size();
    auto& vec = *targets;
    if (greater_than_indices_.empty()) {  // check only less than constraint
      // filter targets
      uint32_t next = 0;
      for (uint32_t i = 0; i < size; ++i) {
        if (vec[i] < less_than) {
          vec[next++] = vec[i];
        }
      }
      vec.resize(next);
    } else {
      VertexID greater_than = 0;  // max of all constraints
      for (auto index : greater_than_indices_) {
        auto val = group.getKeyVal(index);
        greater_than = std::max(greater_than, val);
      }

      // filter targets
      uint32_t next = 0;
      for (uint32_t i = 0; i < size; ++i) {
        if (vec[i] > greater_than && vec[i] < less_than) {
          vec[next++] = vec[i];
        }
      }
      vec.resize(next);
    }
  }

  bool filter(VertexID target, const CompressedSubgraphs& group) const {
    for (auto index : less_than_indices_) {
      auto val = group.getKeyVal(index);
      if (target >= val) return true;
    }
    for (auto index : greater_than_indices_) {
      auto val = group.getKeyVal(index);
      if (target <= val) return true;
    }
    return false;
  }

};  // class TargetFilter

}  // namespace circinus

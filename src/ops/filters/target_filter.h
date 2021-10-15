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

#include "utils/utils.h"

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
      : less_than_indices_(std::move(lt_conditions)), greater_than_indices_(std::move(gt_conditions)) {
    // LOG(INFO) << "lt " << toString(less_than_indices_) << " gt " << toString(greater_than_indices_);
  }

  void filter(std::vector<VertexID>* targets, const CompressedSubgraphs& group) const;

  bool filter(VertexID target, const CompressedSubgraphs& group) const;

  // TODO(tatiana): filter CandidateSetView and return new range?

};  // class TargetFilter

}  // namespace circinus

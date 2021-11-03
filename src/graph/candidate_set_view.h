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

#include <iterator>
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"
#include "graph/vertex_set_view.h"

namespace circinus {

class CandidateSetView : public VertexSetView {
 public:
  CandidateSetView() {}
  using ConstIterator = const VertexID*;

  explicit CandidateSetView(const std::vector<VertexID>& candidates) {
    addRange(candidates.data(), candidates.data() + candidates.size());
  }

  CandidateSetView(const std::vector<VertexID>* candidates, const CandidateScope& scope,
                   const std::vector<VertexID>& partition_offsets);

  CandidateSetView(const std::vector<std::vector<VertexID>>& partitioned_candidates, const CandidateScope& scope);

  inline VertexID operator[](uint32_t index) const {
    DCHECK(!ranges_.empty());
    DCHECK_LT(index, ranges_.front().second);
    return *(ranges_.front().first + index);
  }

  inline ConstIterator begin() const { return ranges_.empty() ? nullptr : ranges_.front().first; }
  inline ConstIterator end() const {
    return ranges_.empty() ? nullptr : (ranges_.front().first + ranges_.front().second);
  }
};

}  // namespace circinus

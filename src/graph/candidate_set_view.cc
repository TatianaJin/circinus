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

#include "graph/candidate_set_view.h"

#include <vector>

#include "glog/logging.h"

#include "graph/types.h"

namespace circinus {

CandidateSetView::CandidateSetView(const std::vector<VertexID>* candidates, const CandidateScope& scope,
                                   const std::vector<VertexID>& partition_offsets) {
  CHECK_EQ(0, partition_offsets.front());
  CHECK_EQ(candidates->size(), partition_offsets.back());
  switch (scope.getType()) {
  case CandidateScopeType::All: {
    if (!candidates->empty()) {
      addRange(candidates->data(), candidates->data() + candidates->size());
    }
    return;
  }
  case CandidateScopeType::Partition: {
    DCHECK_LT(scope.getPartition() + 1, partition_offsets.size());
    addRange(candidates->data() + partition_offsets[scope.getPartition()],
             candidates->data() + partition_offsets[scope.getPartition() + 1]);
    return;
  }
  case CandidateScopeType::Inverse: {
    if (scope.getPartition() == 0) {
      addRange(candidates->data() + partition_offsets[1], candidates->data() + candidates->size());
    } else if (scope.getPartition() + 2 == partition_offsets.size()) {
      addRange(candidates->data(), candidates->data() + partition_offsets[scope.getPartition()]);
    } else {  // if the excluded partition is not at either end, two ranges
      addRange(candidates->data(), candidates->data() + partition_offsets[scope.getPartition()]);
      addRange(candidates->data() + partition_offsets[scope.getPartition() + 1],
               candidates->data() + candidates->size());
    }
  }
  }
}

}  // namespace circinus

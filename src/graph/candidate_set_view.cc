#include "graph/candidate_set_view.h"

#include <vector>

#include "glog/logging.h"

#include "graph/types.h"

namespace circinus {

CandidateSetView::CandidateSetView(const std::vector<VertexID>* candidates, const CandidateScope& scope,
                                   const std::vector<VertexID>& partition_offsets) {
  DCHECK(candidates != nullptr);
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
    return;
  }
  case CandidateScopeType::PartitionRange: {
    addRange(candidates->data() + partition_offsets[scope.getPartition()] + scope.getRangeStart(),
             candidates->data() + partition_offsets[scope.getPartition()] + scope.getRangeEnd() + 1);
    return;
  }
  case CandidateScopeType::Range: {
    addRange(candidates->data() + scope.getRangeStart(), candidates->data() + scope.getRangeEnd() + 1);
  }
  }
  CHECK_EQ(ranges_.size(), 1) << "For CandidateSetView, only 1 range is supported for efficiency reason";
}

CandidateSetView::CandidateSetView(const std::vector<std::vector<VertexID>>& partitioned_candidates,
                                   const CandidateScope& scope) {
  switch (scope.getType()) {
  case CandidateScopeType::All: {
    for (auto& range : partitioned_candidates) {
      addRange(range.data(), range.data() + range.size());
    }
    return;
  }
  case CandidateScopeType::Partition: {
    DCHECK_LT(scope.getPartition(), partitioned_candidates.size());
    addRange(partitioned_candidates[scope.getPartition()].data(),
             partitioned_candidates[scope.getPartition()].data() + partitioned_candidates[scope.getPartition()].size());
    return;
  }
  case CandidateScopeType::Inverse: {
    for (uint32_t partition = 0; partition < partitioned_candidates.size(); ++partition) {
      if (partition == scope.getPartition()) continue;
      addRange(partitioned_candidates[partition].data(),
               partitioned_candidates[partition].data() + partitioned_candidates[partition].size());
    }
    return;
  }
  case CandidateScopeType::PartitionRange: {
    addRange(partitioned_candidates[scope.getPartition()].data() + scope.getRangeStart(),
             partitioned_candidates[scope.getPartition()].data() + scope.getRangeEnd() + 1);
  }
  default: { scope.print(LOG(FATAL) << "Unsupported scope type "); }
  }

  CHECK_EQ(ranges_.size(), 1) << "For CandidateSetView, only 1 range is supported for efficiency reason";
}

}  // namespace circinus

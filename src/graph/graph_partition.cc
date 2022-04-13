#include "graph/graph_partition.h"

#include <memory>

namespace circinus {

std::unique_ptr<GraphPartitionBase> GraphPartitionBase::createGraphPartition(const CandidateScope& src_scope,
                                                                             const CandidateScope& dst_scope,
                                                                             const ReorderedPartitionedGraph* graph) {
  if (src_scope.getType() == CandidateScopeType::Partition) {
    if (dst_scope.getType() == CandidateScopeType::Partition) {
      return std::make_unique<SDGraphPartition<true, true>>(graph, src_scope.getPartition(), dst_scope.getPartition());
    }
    return std::make_unique<SDGraphPartition<true, false>>(graph, src_scope.getPartition());
  }
  if (dst_scope.getType() == CandidateScopeType::Partition) {
    return std::make_unique<SDGraphPartition<false, true>>(graph, 0, dst_scope.getPartition());
  }
  return std::make_unique<SDGraphPartition<false, false>>(graph, 0);
}

}  // namespace circinus

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

#include "graph/partitioned_graph.h"

namespace circinus {

/**
 * A partition view of ReorderedPartitionedGraph.
 */
class GraphPartition {
  const ReorderedPartitionedGraph* const original_graph_ = nullptr;
  const uint32_t src_partition_;  // the major partition
  const uint32_t dst_partition_;  // the (local) partition that shows a view of the neighborhood

 public:
  GraphPartition(const ReorderedPartitionedGraph* original_graph, uint32_t partition)
      : GraphPartition(original_graph, partition, partition) {}

  GraphPartition(const ReorderedPartitionedGraph* original_graph, uint32_t src_partition, uint32_t dst_partition)
      : original_graph_(original_graph), src_partition_(src_partition), dst_partition_(dst_partition) {}

  inline VertexID getPartitionOffset() const { return original_graph_->getPartitionOffset(src_partition_); }
  inline VertexID getPartitionEnd() const { return original_graph_->getPartitionEnd(src_partition_); }

  inline VertexID getOriginalId(VertexID id) const { return original_graph_->getOriginalVertexId(id); }

  /* start of vertex accessors within src-dst partition */

  /** @returns The range of (reordered) vertex ids within src_partition_. */
  inline const std::pair<VertexID, VertexID>& getVertexRangeByLabel(LabelID lid) const {
    return original_graph_->getVertexRangeByLabel(lid, src_partition_);
  }

  /* end of vertex accessors within src-dst partition */

  /* start of full neighborhood accessors for filtering */

  /** @param vid The reordered vertex id. */
  inline LabelID getVertexLabel(VertexID vid) const {
    if (vid >= getPartitionOffset() && vid < getPartitionEnd()) {
      return original_graph_->getVertexLabelInPartition(vid, src_partition_);
    }
    return original_graph_->getVertexLabel(vid);
  }

  /** @param vid The reordered vertex id. */
  inline VertexID getVertexOutDegree(VertexID vid) const { return original_graph_->getVertexOutDegree(vid); }

  /** @param vid The reordered vertex id. */
  inline auto getOutNeighbors(VertexID vid) const { return original_graph_->getOutNeighbors(vid); }

  /* end of full neighborhood accessors for filtering */
};

}  // namespace circinus

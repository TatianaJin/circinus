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

#include <memory>
#include <utility>

#include "graph/partitioned_graph.h"

namespace circinus {

class GraphPartitionBase {
 public:
  using NeighborSet = VertexSetView;

  virtual ~GraphPartitionBase() {}
  virtual bool containsDestination(VertexID vid) const = 0;
  virtual NeighborSet getOutNeighborsWithHint(VertexID vid, LabelID nbr_label) const = 0;
  virtual NeighborSet getInNeighborsWithHint(VertexID vid, LabelID nbr_label) const = 0;
  virtual VertexID getVertexOutDegreeWithHint(VertexID vid, LabelID nbr_label) const = 0;
  virtual VertexID getVertexInDegreeWithHint(VertexID vid, LabelID nbr_label) const = 0;
  static std::unique_ptr<GraphPartitionBase> createGraphPartition(const CandidateScope& src_scope,
                                                                  const CandidateScope& dst_scope,
                                                                  const ReorderedPartitionedGraph* graph);
};

/**
 * A partition view of ReorderedPartitionedGraph.
 *
 * The following functions are implemented for backtracking
 *
 * Check if the vertex is in dst partition(s)
 * bool containsDestination(VertexID vid) const;
 *
 * Get neighbors in dst partition(s):
 * VertexSetView getOutNeighborsWithHint(VertexID vid, LabelID nbr_label) const;
 *
 * Get neighbors in src partition(s):
 * VertexSetView getInNeighborsWithHint(VertexID vid, LabelID nbr_label) const;
 *
 * Get the number (upper bound) of neighbors in dst partition(s):
 * VertexID getVertexOutDegreeWithHint(VertexID vid, LabelID nbr_label) const;
 *
 * Get the number (upper bound) of neighbors in src partition(s):
 * VertexID getVertexInDegreeWithHint(VertexID vid, LabelID nbr_label) const;
 */
template <bool limit_src, bool limit_dst>
class SDGraphPartition : public GraphPartitionBase {
 protected:
  const ReorderedPartitionedGraph* const original_graph_ = nullptr;
  const uint32_t src_partition_;
  const uint32_t dst_partition_;

 public:
  using NeighborSet = VertexSetView;

  SDGraphPartition(const ReorderedPartitionedGraph* original_graph, uint32_t src_partition, uint32_t dst_partition)
      : original_graph_(original_graph), src_partition_(src_partition), dst_partition_(dst_partition) {}

  SDGraphPartition(const ReorderedPartitionedGraph* original_graph, uint32_t src_partition)
      : SDGraphPartition(original_graph, src_partition, 0) {}

  inline VertexID getPartitionOffset() const { return original_graph_->getPartitionOffset(src_partition_); }
  inline VertexID getPartitionEnd() const { return original_graph_->getPartitionEnd(src_partition_); }

  inline VertexID getOriginalId(VertexID id) const { return original_graph_->getOriginalVertexId(id); }

  /* start of vertex accessors within src-dst partition */

  /** @returns The range of (reordered) vertex ids within src_partition_. */
  inline std::pair<VertexID, VertexID> getVertexRangeByLabel(LabelID lid) const {
    return original_graph_->getVertexRangeByLabel(lid, src_partition_);
  }

  inline bool containsDestination(VertexID vid) const override {
    if
      constexpr(limit_dst) {
        return vid >= original_graph_->getPartitionOffset(dst_partition_) &&
               vid < original_graph_->getPartitionEnd(dst_partition_);
      }
    return vid < original_graph_->getNumVertices();
  }

  inline NeighborSet getOutNeighborsWithHint(VertexID vid, LabelID nbr_label) const override {
    if
      constexpr(limit_dst) {
        return original_graph_->getOutNeighborsInPartitionWithHint(vid, nbr_label, dst_partition_);
      }
    return original_graph_->getAllOutNeighborsWithHint(vid, nbr_label);
  }

  inline NeighborSet getInNeighborsWithHint(VertexID vid, LabelID nbr_label) const override {
    if
      constexpr(limit_src) {
        return original_graph_->getOutNeighborsInPartitionWithHint(vid, nbr_label, src_partition_);
      }
    return original_graph_->getAllOutNeighborsWithHint(vid, nbr_label);
  }

  inline VertexID getVertexOutDegreeWithHint(VertexID vid, LabelID nbr_label) const override {
    if
      constexpr(limit_dst) {
        return original_graph_->getVertexOutDegreeInPartitionWithHint(vid, nbr_label, dst_partition_);
      }
    return original_graph_->getVertexOutDegree(vid);
  }

  inline VertexID getVertexInDegreeWithHint(VertexID vid, LabelID nbr_label) const override {
    if
      constexpr(limit_src) {
        return original_graph_->getVertexOutDegreeInPartitionWithHint(vid, nbr_label, src_partition_);
      }
    return original_graph_->getVertexOutDegree(vid);
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

using GraphPartition = SDGraphPartition<true, false>;

}  // namespace circinus

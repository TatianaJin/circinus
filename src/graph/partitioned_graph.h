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

#include <algorithm>
#include <numeric>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "graph/graph.h"
#include "graph/graph_base.h"
#include "graph/types.h"
#include "graph/vertex_set_view.h"

namespace circinus {

/* Reordered and partitioned graph
 *
 * The vertices are first ordered by their partitions, then by their labels, and then optionally ordered by their
 * degrees in descending order.
 */
class ReorderedPartitionedGraph : public GraphBase {
  static const std::pair<VertexID, VertexID> ZERO_RANGE;

  /* statistics */
  uint32_t n_partitions_;
  uint64_t n_edge_cuts_ = 0;

  /* index */
  std::vector<VertexID> vertex_ids_;         // the original vertex id of the i-th vertex in vlist_
  std::vector<VertexID> partition_offsets_;  // size n_partitions_ + 1
  std::vector<unordered_map<LabelID, std::pair<VertexID, VertexID>>>
      label_ranges_per_part_;  // {partition, {label, {start, end}}}

 public:
  explicit ReorderedPartitionedGraph(const std::string& path, uint32_t n_partitions = 20, bool sort_by_degree = true);

  explicit ReorderedPartitionedGraph(const Graph& graph, uint32_t n_partitions = 20, bool sort_by_degree = true);

  ReorderedPartitionedGraph() {}  // for loading from binary

  inline uint32_t getNumPartitions() const { return n_partitions_; }
  inline uint64_t getNumEdgeCuts() const { return n_edge_cuts_; }

  inline VertexID getOriginalVertexId(VertexID id) const { return vertex_ids_[id]; }

  /* compatible interface for graph partition TODO(tatiana): unnecessary? */
  inline VertexID getVertexGlobalId(VertexID id) const { return id; }
  inline VertexID getVertexLocalId(VertexID id) const { return id; }

  bool containsDestination(VertexID vid) const { return vid < getNumVertices(); }

  /** Get the number of neighbors refined by hints.
   * @param partition The desired partition of neighbors.
   */
  VertexID getVertexOutDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t partition = 0) const {
    auto nbrs = getOutNeighbors(id);
    auto first_neighbor = nbrs.first[0];
    auto last_neighbor = nbrs.first[nbrs.second - 1];
    // get the range of vertex ids that satisfy the hint
    VertexID range_l, range_r;
    if (nbr_label == ALL_LABEL) {
      range_l = partition_offsets_[partition];
      range_r = partition_offsets_[partition + 1];
    } else {
      std::tie(range_l, range_r) = label_ranges_per_part_[partition].at(nbr_label);
    }
    if (first_neighbor >= range_r || last_neighbor < range_l) {
      return 0;
    }
    // upper bound is (range_r - range_l)
    range_r = std::min(range_r, last_neighbor + 1);
    range_l = std::max(range_l, first_neighbor);
    return std::min(range_r - range_l, getVertexOutDegree(id));
  }

  inline VertexID getVertexInDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t partition = 0) const {
    return getVertexOutDegreeWithHint(id, nbr_label, partition);
  }

  // do we need to consider an O(1) time implementation?
  inline LabelID getVertexLabel(VertexID id) const { return getVertexLabelInPartition(id, getVertexPartition(id)); }

  LabelID getVertexLabelInPartition(VertexID id, uint32_t partition) const {
    auto& label_ranges = label_ranges_per_part_[partition];
    for (auto& pair : label_ranges) {
      if (pair.second.first <= id && id < pair.second.second) {
        return pair.first;
      }
    }
    LOG(FATAL) << "cannot find label range for vertex " << id << " partition " << partition;
    return 0;
  }

  static std::pair<const VertexID*, const VertexID*> getVertexRange(const VertexID* search_start,
                                                                    const VertexID* search_end, VertexID range_l,
                                                                    VertexID range_r) {
    const VertexID* start = nullptr;
    const VertexID* end = nullptr;
    if (std::distance(search_start, search_end) >= 32) {  // binary search
      start = std::lower_bound(search_start, search_end, range_l);
      if (start < search_end) {
        end = std::lower_bound(start, search_end, range_r);
        return {start, end};
      }
      return {search_end, search_end};
    }
    // linear scan
    end = start = search_end;
    for (auto ptr = search_start; ptr < search_end; ++ptr) {
      if (*ptr >= range_l) {
        start = ptr;
      }
      if (*ptr >= range_r) {
        end = ptr;
        break;
      }
    }
    return {start, end};
  }

  /** Get the neighbors that satisfy the hints.
   * @param partition The desired partition of neighbors.
   */
  VertexSetView getOutNeighborsWithHint(VertexID id, LabelID nbr_label, uint32_t partition = 0) const {
    auto nbrs = getOutNeighbors(id);
    VertexID range_l, range_r;
    // get the range of vertex ids that satisfy the hint
    if (nbr_label == ALL_LABEL) {
      range_l = partition_offsets_[partition];
      range_r = partition_offsets_[partition + 1];
    } else {
      std::tie(range_l, range_r) = label_ranges_per_part_[partition].at(nbr_label);
    }
    auto[start, end] = getVertexRange(nbrs.first, nbrs.first + nbrs.second, range_l, range_r);
    return VertexSetView(start, end);
  }

  inline VertexSetView getInNeighborsWithHint(VertexID id, LabelID nbr_label, uint32_t partition = 0) const {
    return getOutNeighborsWithHint(id, nbr_label, partition);
  }

  VertexSetView getAllOutNeighborsWithHint(VertexID id, LabelID nbr_label) const {
    VertexSetView view;
    auto nbrs = getOutNeighbors(id);
    if (nbr_label == ALL_LABEL) {
      view.addRange(nbrs.first, nbrs.first + nbrs.second);
      return view;
    }
    const VertexID* search_start = nbrs.first;
    const VertexID* const search_end = nbrs.first + nbrs.second;
    for (uint32_t partition = 0; partition < getNumPartitions(); ++partition) {
      // get the range of vertex ids that satisfy the hint
      auto[range_l, range_r] = label_ranges_per_part_[partition].at(nbr_label);
      if (range_l == range_r) continue;
      auto[start, end] = getVertexRange(search_start, search_end, range_l, range_r);
      if (start == search_end) {
        break;
      }
      view.addRange(start, end);
      search_start = end;
    }
    return view;
  }

  uint32_t getVertexPartition(VertexID id) const {
    uint32_t partition = 0;
    for (; partition < n_partitions_; ++partition) {
      if (partition_offsets_[partition] > id) {
        break;
      }
    }
    return --partition;
  }

  void clear() override {
    GraphBase::clear();
    vertex_ids_.clear();
    partition_offsets_.clear();
    label_ranges_per_part_.clear();
    n_partitions_ = 0;
    n_edge_cuts_ = 0;
  }

  /* end of overload functions */

  /**
   * @returns The range [first, last) of ids of the vertices in the partition.
   */
  inline std::pair<VertexID, VertexID> getPartitionRange(uint32_t partition) const {
    DCHECK_LT(partition, n_partitions_);
    return std::make_pair(partition_offsets_[partition], partition_offsets_[partition + 1]);
  }

  inline VertexID getPartitionOffset(uint32_t partition) const { return partition_offsets_[partition]; }
  inline VertexID getPartitionEnd(uint32_t partition) const { return partition_offsets_[partition + 1]; }

  inline const auto& getLabelOffsetsInPartition(uint32_t partition) const {
    DCHECK_LT(partition, n_partitions_);
    return label_ranges_per_part_[partition];
  }

  inline const std::pair<VertexID, VertexID>& getVertexRangeByLabel(LabelID lid, uint32_t partition) const {
    auto& label_ranges = label_ranges_per_part_[partition];
    auto pos = label_ranges.find(lid);
    if (pos == label_ranges.end()) {
      return ZERO_RANGE;
    }
    return pos->second;
  }

  /* persistence */
  void dumpToFile(const std::string& path) const override;
  void loadUndirectedGraphFromBinary(std::istream& input) override;

  std::pair<double, double> getMemoryUsage() const override;

 protected:
  void saveAsBinaryInner(std::ostream& output) const override;

  template <bool by_partition, bool by_degree>
  void sortVertices(const std::vector<idx_t>& parts, const std::vector<LabelID> labels, const GraphBase* src_graph);

  void reorder(bool sort_by_degree, const std::vector<LabelID>& labels, const GraphBase*);
  void reorder(bool sort_by_degree, const std::vector<LabelID>& labels, const std::vector<idx_t>& parts,
               const GraphBase*);
  void computeLabelOffsets(const std::vector<LabelID>& labels);
  void reconstructByOrder();
  void reconstructByOrder(const GraphBase& src_graph);
};

}  // namespace circinus

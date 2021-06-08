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
#include <unordered_map>
#include <utility>
#include <vector>

#include "graph/graph.h"
#include "graph/graph_base.h"
#include "graph/types.h"

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

  inline VertexID getOriginalVertexId(VertexID id) const { return vertex_ids_[id]; }

  /* compatible interface for graph partition TODO(tatiana): unnecessary? */
  inline VertexID getVertexGlobalId(VertexID id) const { return id; }
  inline VertexID getVertexLocalId(VertexID id) const { return id; }

  /** Get the number of neighbors refined by hints.
   * @param partition The desired partition of neighbors.
   */
  inline VertexID getVertexOutDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t partition = 0) const {
    auto nbrs = getOutNeighbors(id);
    auto first_neighbor = nbrs.first[0];
    auto last_neighbor = nbrs.first[nbrs.second - 1];
    // get the range of vertex ids that satisfy the hint
    auto[range_l, range_r] = label_ranges_per_part_[partition].at(nbr_label);
    if (first_neighbor >= range_r || last_neighbor < range_l) {
      return 0;
    }
    // upper bound is (range_r - range_l)
    range_r = std::min(range_r, last_neighbor + 1);
    range_l = std::max(range_l, first_neighbor);
    return std::min(range_r - range_l, getVertexOutDegree(id));
  }

  // do we need to consider an O(1) time implementation?
  inline LabelID getVertexLabel(VertexID id) const {
    uint32_t partition = 0;
    for (; partition < n_partitions_; ++partition) {
      if (partition_offsets_[partition] > id) {
        break;
      }
    }
    return getVertexLabelInPartition(id, --partition);
  }

  LabelID getVertexLabelInPartition(VertexID id, uint32_t partition) const {
    auto& label_ranges = label_ranges_per_part_[partition];
    for (auto& pair : label_ranges) {
      if (pair.second.first <= id && id < pair.second.second) {
        return pair.first;
      }
    }
    LOG(ERROR) << "cannot find label range " << id << " partition " << partition;
    return 0;
  }

  /** Get the neighbors that satisfy the hints.
   * @param partition The desired partition of neighbors.
   */
  inline std::pair<const VertexID*, uint32_t> getOutNeighborsWithHint(VertexID id, LabelID nbr_label,
                                                                      uint32_t partition = 0) const {
    auto nbrs = getOutNeighbors(id);
    const VertexID* start = nullptr;
    const VertexID* end = nullptr;
    // get the range of vertex ids that satisfy the hint
    auto[range_l, range_r] = label_ranges_per_part_[partition].at(nbr_label);
    if (getVertexOutDegree(id) >= 32) {  // binary search
      start = std::lower_bound(nbrs.first, nbrs.first + nbrs.second, range_l);
      if (start - nbrs.first < nbrs.second) {
        end = std::lower_bound(start, nbrs.first + nbrs.second, range_r);
        return {start, end - start};
      }
      return {nullptr, 0};
    } else {  // linear scan
      for (uint32_t i = 0; i < nbrs.second; ++i) {
        if (nbrs.first[i] >= range_l) {
          start = nbrs.first + i;
        }
        if (nbrs.first[i] >= range_r) {
          end = nbrs.first + i;
          break;
        }
      }
      return {start, end - start};
    }
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

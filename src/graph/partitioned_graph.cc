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

#include "graph/partitioned_graph.h"

#include <algorithm>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"
#include "utils/file_utils.h"
#include "utils/utils.h"

namespace circinus {

const std::pair<VertexID, VertexID> ReorderedPartitionedGraph::ZERO_RANGE = {0, 0};

ReorderedPartitionedGraph::ReorderedPartitionedGraph(const std::string& path, uint32_t n_partitions,
                                                     bool sort_by_degree)
    : n_partitions_(n_partitions), partition_offsets_(n_partitions + 1, 0), label_offsets_per_part_(n_partitions) {
  auto labels = loadUndirectedGraph(path);
  vertex_ids_.resize(getNumVertices());
  std::iota(vertex_ids_.begin(), vertex_ids_.end(), 0);
  // TODO(tatiana): now only single-thread processing, support parallel processing later
  // partition and reorder
  if (n_partitions > 1) {
    auto[parts, n_edge_cuts] = getMetisParts(n_partitions);
    n_edge_cuts_ = n_edge_cuts;
    reorder(sort_by_degree, labels, parts, this);
  } else {
    reorder(sort_by_degree, labels, this);
  }
  computeLabelOffsets(labels);
  reconstructByOrder();
}

ReorderedPartitionedGraph::ReorderedPartitionedGraph(const std::string& path, const std::string& part_path,
                                                     uint32_t n_partitions, bool sort_by_degree)
    : n_partitions_(n_partitions), partition_offsets_(n_partitions + 1, 0), label_offsets_per_part_(n_partitions) {
  auto labels = loadUndirectedGraph(path);
  vertex_ids_.resize(getNumVertices());
  std::iota(vertex_ids_.begin(), vertex_ids_.end(), 0);
  // read partition and reorder
  auto parts = getFennelParts(part_path);
  reorder(sort_by_degree, labels, parts, this);
  computeLabelOffsets(labels);
  reconstructByOrder();
}

ReorderedPartitionedGraph::ReorderedPartitionedGraph(const Graph& graph, uint32_t n_partitions, bool sort_by_degree)
    : n_partitions_(n_partitions),
      vertex_ids_(graph.getNumVertices()),
      partition_offsets_(n_partitions + 1, 0),
      label_offsets_per_part_(n_partitions) {
  copyMetadata(graph);
  auto& labels = graph.getVertexLabels();
  std::iota(vertex_ids_.begin(), vertex_ids_.end(), 0);
  if (n_partitions > 1) {
    auto[parts, n_edge_cuts] = graph.getMetisParts(n_partitions);
    n_edge_cuts_ = n_edge_cuts;
    reorder(sort_by_degree, labels, parts, &graph);
  } else {
    reorder(sort_by_degree, labels, &graph);
  }
  computeLabelOffsets(labels);
  reconstructByOrder(graph);
}

template <bool by_partition, bool by_degree>
void ReorderedPartitionedGraph::sortVertices(const std::vector<idx_t>& parts, const std::vector<LabelID> labels,
                                             const GraphBase* src_graph) {
  std::sort(vertex_ids_.begin(), vertex_ids_.end(), [src_graph, &labels, &parts](VertexID v1, VertexID v2) {
    if (!by_partition || parts[v1] == parts[v2]) {
      if (labels[v1] == labels[v2]) {
        return (!by_degree && v1 < v2) || src_graph->getVertexOutDegree(v1) < src_graph->getVertexOutDegree(v2);
      }
      return labels[v1] < labels[v2];
    }
    return parts[v1] < parts[v2];
  });
}

// TODO(tatiana): now only single-thread processing, support parallel processing later
void ReorderedPartitionedGraph::reorder(bool sort_by_degree, const std::vector<LabelID>& labels,
                                        const std::vector<idx_t>& parts, const GraphBase* src_graph) {
  if (sort_by_degree) {
    sortVertices<true, true>(parts, labels, src_graph);
  } else {
    sortVertices<true, false>(parts, labels, src_graph);
  }

  // construct partition offsets
  auto start = vertex_ids_.begin();
  for (uint32_t i = 1; i < n_partitions_; ++i) {
    start = std::lower_bound(start, vertex_ids_.end(), i,
                             [this, &parts](VertexID v, uint32_t partition) { return parts[v] < partition; });
    DCHECK_EQ(parts[*start], i);
    partition_offsets_[i] = start - vertex_ids_.begin();
  }
  partition_offsets_.back() = vertex_ids_.size();
}

void ReorderedPartitionedGraph::reorder(bool sort_by_degree, const std::vector<LabelID>& labels,
                                        const GraphBase* src_graph) {
  if (sort_by_degree) {
    sortVertices<false, true>({}, labels, src_graph);
  } else {
    sortVertices<false, false>({}, labels, src_graph);
  }
  partition_offsets_ = {0, getNumVertices()};
}

void ReorderedPartitionedGraph::computeLabelOffsets(const std::vector<LabelID>& labels) {
  std::set<LabelID> label_set;
  for (auto& pair : vertex_cardinality_by_label_) {
    label_set.insert(pair.first);
  }
  LabelID n_labels = label_set.size();
  label_idx_.reserve(n_labels);
  label_set_.resize(n_labels);
  LabelID idx = 0;
  for (auto l : label_set) {
    label_set_[idx] = l;
    label_idx_.insert({l, idx++});
  }

  auto start = vertex_ids_.begin();
  for (uint32_t i = 0; i < n_partitions_; ++i) {
    label_offsets_per_part_[i].resize(label_idx_.size() + 1, INVALID_VERTEX_ID);
    auto end = vertex_ids_.begin() + partition_offsets_[i + 1];

    label_offsets_per_part_[i][0] = partition_offsets_[i];
    LabelID last_label_idx = 0;  // last label whose offset is set
    VertexID last_offset = partition_offsets_[i];
    while (start < end) {
      LabelID current_label = labels[*start];
      auto range_end = std::lower_bound(start + 1, end, current_label + 1,
                                        [&labels](VertexID v, LabelID label) { return labels[v] < label; });
      auto label_idx = label_idx_.at(labels[*start]);
      while (last_label_idx < label_idx) {
        label_offsets_per_part_[i][++last_label_idx] = last_offset;
      }
      last_label_idx = label_idx + 1;
      last_offset = label_offsets_per_part_[i][last_label_idx] = range_end - vertex_ids_.begin();
      start = range_end;
    }
    while (last_label_idx < label_idx_.size()) {
      label_offsets_per_part_[i][++last_label_idx] = last_offset;
    }
  }
}

void ReorderedPartitionedGraph::reconstructByOrder() {
  // reconstruct vlist and elist
  std::vector<VertexID> vertex_order(
      vertex_ids_.size());  // original vertex id, the index in vlist_ of the vertex with the original id
  for (VertexID i = 0; i < vertex_ids_.size(); ++i) {
    vertex_order[vertex_ids_[i]] = i;
  }
  std::vector<VertexID> new_elist(elist_.size());
  std::vector<EdgeID> new_vlist(vlist_.size());
  VertexID new_vertex_id = 0;
  new_vlist[new_vertex_id] = 0;
  for (auto v : vertex_ids_) {
    new_vlist[new_vertex_id + 1] = new_vlist[new_vertex_id];
    auto nbrs = getOutNeighbors(v);
    for (uint32_t i = 0; i < nbrs.second; ++i) {
      new_elist[new_vlist[new_vertex_id + 1]] = vertex_order[nbrs.first[i]];
      ++new_vlist[new_vertex_id + 1];
    }
    ++new_vertex_id;
  }
  vlist_.swap(new_vlist);
  elist_.swap(new_elist);
  // sort neighbors by id for each vertex
  for (VertexID i = 0; i < n_vertices_; ++i) {
    std::sort(elist_.begin() + vlist_[i], elist_.begin() + vlist_[i + 1]);
  }
}

void ReorderedPartitionedGraph::reconstructByOrder(const GraphBase& src_graph) {
  // reconstruct vlist and elist
  std::vector<VertexID> vertex_order(
      vertex_ids_.size());  // original vertex id, the index in vlist_ of the vertex with the original id
  for (VertexID i = 0; i < vertex_ids_.size(); ++i) {
    vertex_order[vertex_ids_[i]] = i;
  }
  elist_.resize(src_graph.getNumEdges() * 2);
  vlist_.resize(src_graph.getNumVertices() + 1);
  VertexID new_vertex_id = 0;
  vlist_[new_vertex_id] = 0;
  for (auto v : vertex_ids_) {
    vlist_[new_vertex_id + 1] = vlist_[new_vertex_id];
    auto nbrs = src_graph.getOutNeighbors(v);
    for (uint32_t i = 0; i < nbrs.second; ++i) {
      elist_[vlist_[new_vertex_id + 1]] = vertex_order[nbrs.first[i]];
      ++vlist_[new_vertex_id + 1];
    }
    ++new_vertex_id;
  }
  // sort neighbors by id for each vertex
  for (VertexID i = 0; i < n_vertices_; ++i) {
    std::sort(elist_.begin() + vlist_[i], elist_.begin() + vlist_[i + 1]);
  }
}

void ReorderedPartitionedGraph::dumpToFile(const std::string& path) const {
  std::vector<VertexID> vertex_order(
      vertex_ids_.size());  // original vertex id, the index in vlist_ of the vertex with the original id
  for (VertexID i = 0; i < vertex_ids_.size(); ++i) {
    vertex_order[vertex_ids_[i]] = i;
  }

  auto output = openOutputFile(path);
  output << "t " << getNumVertices() << ' ' << getNumEdges() << '\n';

  for (VertexID i = 0; i < getNumVertices(); ++i) {
    output << "v " << i << ' ' << getVertexLabel(vertex_order[i]) << ' ' << getVertexOutDegree(vertex_order[i]) << '\n';
  }

  for (VertexID i = 0; i < getNumVertices(); ++i) {
    auto neighbors = getOutNeighbors(vertex_order[i]);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (i < vertex_ids_[neighbors.first[j]]) {
        output << "e " << i << ' ' << vertex_ids_[neighbors.first[j]] << '\n';
      }
    }
  }
  output.close();
}

void ReorderedPartitionedGraph::saveAsBinaryInner(std::ostream& output) const {
  bool partitioned_graph = true;
  output.write(reinterpret_cast<const char*>(&partitioned_graph), sizeof(bool));
  // put at the beginning to faliciate parallel loading when needed
  output.write(reinterpret_cast<const char*>(&n_partitions_), sizeof(n_partitions_));

  GraphBase::saveAsBinaryInner(output);
  output.write(reinterpret_cast<const char*>(&n_edge_cuts_), sizeof(n_edge_cuts_));

  vectorToBinaryStream(output, vertex_ids_);
  vectorToBinaryStream(output, partition_offsets_);

  {  // label_idx_
    auto size = label_idx_.size();
    output.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& pair : label_idx_) {
      output.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
      output.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second));
    }
  }
  vectorToBinaryStream(output, label_set_);  // label_set_
  {                                          // label ranges per partition
    auto size = label_offsets_per_part_.size();
    output.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& vec : label_offsets_per_part_) {
      vectorToBinaryStream(output, vec);
    }
  }
}

void ReorderedPartitionedGraph::loadUndirectedGraphFromBinary(std::istream& input) {
  uint32_t n_partitions;
  input.read(reinterpret_cast<char*>(&n_partitions), sizeof(n_partitions));
  GraphBase::loadUndirectedGraphFromBinary(input);
  n_partitions_ = n_partitions;
  input.read(reinterpret_cast<char*>(&n_edge_cuts_), sizeof(n_edge_cuts_));
  LOG(INFO) << "n_edge_cuts_ = " << n_edge_cuts_;

  binaryStreamToVector(input, vertex_ids_);
  CHECK_EQ(vertex_ids_.size(), getNumVertices());
  binaryStreamToVector(input, partition_offsets_);
  CHECK_EQ(partition_offsets_.size(), n_partitions_ + 1);

  {  // label_idx_
    size_t size;
    input.read(reinterpret_cast<char*>(&size), sizeof(size));
    label_idx_.reserve(size);
    for (size_t i = 0; i < size; ++i) {
      LabelID label, idx;
      input.read(reinterpret_cast<char*>(&label), sizeof(LabelID));
      input.read(reinterpret_cast<char*>(&idx), sizeof(LabelID));
      label_idx_.insert({label, idx});
    }
  }
  binaryStreamToVector(input, label_set_);  // label_set_
  {                                         // label ranges per partition
    size_t size;
    input.read(reinterpret_cast<char*>(&size), sizeof(size));
    label_offsets_per_part_.resize(size);
    for (size_t j = 0; j < size; ++j) {
      binaryStreamToVector(input, label_offsets_per_part_[j]);
    }
  }
}

std::pair<double, double> ReorderedPartitionedGraph::getMemoryUsage() const {
  auto ret = GraphBase::getMemoryUsage();
  ret.first += sizeof(ReorderedPartitionedGraph) - sizeof(GraphBase);
  ret.second += sizeof(ReorderedPartitionedGraph) - sizeof(GraphBase);
  ret.first += vertex_ids_.capacity() * sizeof(VertexID);
  ret.second += vertex_ids_.size() * sizeof(VertexID);
  ret.first += partition_offsets_.capacity() * sizeof(VertexID);
  ret.second += partition_offsets_.size() * sizeof(VertexID);

  auto pair_size = sizeof(unordered_map<LabelID, LabelID>::value_type);
  ret.first += label_idx_.bucket_count() * (pair_size + 1);
  ret.second += label_idx_.bucket_count() * (pair_size + 1);
  ret.first += label_set_.capacity() * sizeof(LabelID);
  ret.second += label_set_.size() * sizeof(LabelID);

  auto vec_size = sizeof(std::vector<VertexID>);
  ret.first += label_offsets_per_part_.capacity() * vec_size;
  ret.second += label_offsets_per_part_.size() * vec_size;
  for (auto& vec : label_offsets_per_part_) {
    ret.first += vec.capacity() * sizeof(VertexID);
    ret.second += vec.size() * sizeof(VertexID);
  }
  return ret;
}

std::vector<EdgeID> ReorderedPartitionedGraph::computeEdgeCounts(uint32_t partition) const {
  std::vector<EdgeID> edge_counts(n_partitions_ - partition, 0);
  EdgeID& intra_edge_count = edge_counts[0];
  for (auto vid = partition_offsets_[partition]; vid < partition_offsets_[partition + 1]; ++vid) {
    auto nbrs = getOutNeighbors(vid);
    for (uint32_t i = 0; i < nbrs.second; ++i) {
      auto nbr = nbrs.first[i];
      if (nbr > vid) {
        if (nbr < partition_offsets_[partition + 1]) {  // intra-partition edge
          ++intra_edge_count;
        } else {  // inter-partition edge
          for (auto nbr_partition = partition + 1; nbr_partition < n_partitions_; ++nbr_partition) {
            if (nbr < partition_offsets_[nbr_partition + 1]) {
              ++edge_counts[nbr_partition - partition];
              break;
            }
          }
        }
      }
    }
  }
  return edge_counts;
}

void ReorderedPartitionedGraph::showPartitionInfo() const {
  LOG(INFO) << "n_partitions\t" << n_partitions_;
  {  // vertex balance
    VertexID min = n_vertices_, max = 0, square_sum = 0;
    for (uint32_t partition = 0; partition < n_partitions_; ++partition) {
      auto partition_size = partition_offsets_[partition + 1] - partition_offsets_[partition];
      LOG(INFO) << "partition " << partition << " size " << partition_size;
      min = std::min(min, partition_size);
      max = std::max(max, partition_size);
      square_sum += partition_size * partition_size;
    }
    double average = ((double)n_vertices_) / n_partitions_;
    LOG(INFO) << "[vertex count in partitions]\nmin\t" << min << "\nmax\t" << max << "\navg\t" << average << "\nvar\t"
              << (((double)square_sum) / n_partitions_ - average * average);
  }
  {  // edge balance
    std::vector<std::vector<EdgeID>> edge_counts(n_partitions_);
    std::vector<std::thread> threads;
    for (uint32_t partition = 0; partition < n_partitions_; ++partition) {
      threads.emplace_back(
          [&edge_counts, partition, this]() { edge_counts[partition] = computeEdgeCounts(partition); });
    }
    for (auto& t : threads) t.join();
    EdgeID min = n_edges_, max = 0, square_sum = 0, sum = 0;
    for (uint32_t partition = 0; partition < n_partitions_; ++partition) {
      min = std::min(min, edge_counts[partition][0]);
      max = std::max(max, edge_counts[partition][0]);
      square_sum += edge_counts[partition][0] * edge_counts[partition][0];
      sum += edge_counts[partition][0];
    }
    double average = ((double)sum) / n_partitions_;
    LOG(INFO) << "[intra-partition edge count]\nmin\t" << min << "\nmax\t" << max << "\navg\t" << average << "\nvar\t"
              << (((double)square_sum) / n_partitions_ - average * average);
    std::vector<uint32_t> partition_connectivity(n_partitions_);
    for (uint32_t partition = 0; partition < n_partitions_; ++partition) {
      for (uint32_t i = 1; i < edge_counts[partition].size(); ++i) {
        if (edge_counts[partition][i] > 0) {
          ++partition_connectivity[partition];
          ++partition_connectivity[partition + i];
          // LOG(INFO) << partition << " -(" << edge_counts[partition][i] << ")-> " << (partition + i);
        }
      }
    }
    for (uint32_t partition = 0; partition < n_partitions_; ++partition) {
      LOG(INFO) << "partition " << partition << " connectivity " << partition_connectivity[partition];
    }
  }
}

}  // namespace circinus

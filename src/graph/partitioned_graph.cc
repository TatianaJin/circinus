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
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"
#include "utils/file_utils.h"

namespace circinus {

const std::pair<VertexID, VertexID> ReorderedPartitionedGraph::ZERO_RANGE = {0, 0};

template <bool by_partition, bool by_degree>
void ReorderedPartitionedGraph::sortVertices(const std::vector<idx_t>& parts, const std::vector<LabelID> labels) {
  std::sort(vertex_ids_.begin(), vertex_ids_.end(), [this, &labels, &parts](VertexID v1, VertexID v2) {
    if (!by_partition || parts[v1] == parts[v2]) {
      if (labels[v1] == labels[v2]) {
        return !by_degree || getVertexOutDegree(v1) <= getVertexOutDegree(v2);
      }
      return labels[v1] < labels[v2];
    }
    return parts[v1] < parts[v2];
  });
}

// TODO(tatiana): now only single-thread processing, support parallel processing later
void ReorderedPartitionedGraph::reorder(bool sort_by_degree, const std::vector<LabelID>& labels,
                                        const std::vector<idx_t>& parts) {
  if (sort_by_degree) {
    sortVertices<true, true>(parts, labels);
  } else {
    sortVertices<true, false>(parts, labels);
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

void ReorderedPartitionedGraph::reorder(bool sort_by_degree, const std::vector<LabelID>& labels) {
  if (sort_by_degree) {
    sortVertices<false, true>({}, labels);
  } else {
    sortVertices<false, false>({}, labels);
  }
  partition_offsets_ = {0, getNumVertices()};
}

void ReorderedPartitionedGraph::computeLabelOffsets(const std::vector<LabelID>& labels) {
  auto start = vertex_ids_.begin();
  for (uint32_t i = 0; i < n_partitions_; ++i) {
    auto end = vertex_ids_.begin() + partition_offsets_[i + 1];

    while (start < end) {
      LabelID current_label = labels[*start];
      auto range_end = std::lower_bound(start + 1, end, current_label + 1,
                                        [this, &labels](VertexID v, LabelID label) { return labels[v] < label; });
      label_ranges_per_part_[i].insert(
          {labels[*start], {start - vertex_ids_.begin(), range_end - vertex_ids_.begin()}});
      start = range_end;
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
  GraphBase::saveAsBinaryInner(output);
  output.write(reinterpret_cast<const char*>(&n_partitions_), sizeof(n_partitions_));
  output.write(reinterpret_cast<const char*>(&n_edge_cuts_), sizeof(n_edge_cuts_));

  vectorToBinaryStream(output, vertex_ids_);
  vectorToBinaryStream(output, partition_offsets_);

  {  // label ranges per partition
    auto size = label_ranges_per_part_.size();
    output.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& map : label_ranges_per_part_) {
      auto size = map.size();
      output.write(reinterpret_cast<const char*>(&size), sizeof(size));
      for (const auto& pair : map) {
        output.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
        output.write(reinterpret_cast<const char*>(&pair.second.first), sizeof(pair.second.first));
        output.write(reinterpret_cast<const char*>(&pair.second.second), sizeof(pair.second.second));
      }
    }
  }
}

void ReorderedPartitionedGraph::loadUndirectedGraphFromBinary(std::istream& input) {
  GraphBase::loadUndirectedGraphFromBinary(input);
  input.read(reinterpret_cast<char*>(&n_partitions_), sizeof(n_partitions_));
  input.read(reinterpret_cast<char*>(&n_edge_cuts_), sizeof(n_edge_cuts_));

  binaryStreamToVector(input, vertex_ids_);
  CHECK_EQ(vertex_ids_.size(), getNumVertices());
  binaryStreamToVector(input, partition_offsets_);
  CHECK_EQ(partition_offsets_.size(), n_partitions_ + 1);

  {  // label ranges per partition
    size_t size;
    input.read(reinterpret_cast<char*>(&size), sizeof(size));
    label_ranges_per_part_.resize(size);
    for (size_t j = 0; j < size; ++j) {
      size_t map_size;
      input.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
      auto& map = label_ranges_per_part_[j];
      for (size_t i = 0; i < map_size; ++i) {
        LabelID l;
        VertexID start, end;
        input.read(reinterpret_cast<char*>(&l), sizeof(l));
        input.read(reinterpret_cast<char*>(&start), sizeof(VertexID));
        input.read(reinterpret_cast<char*>(&end), sizeof(VertexID));
        map.insert({l, {start, end}});
      }
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
  auto hashmap_size = sizeof(unordered_map<LabelID, std::pair<VertexID, VertexID>>);
  ret.first += label_ranges_per_part_.capacity() * hashmap_size;
  ret.second += label_ranges_per_part_.size() * hashmap_size;
  auto pair_size = sizeof(unordered_map<LabelID, std::pair<VertexID, VertexID>>::value_type);
  for (auto& hashmap : label_ranges_per_part_) {
    ret.first += hashmap.bucket_count() * (pair_size + 1);
    ret.second += hashmap.bucket_count() * (pair_size + 1);
  }
  return ret;
}

}  // namespace circinus

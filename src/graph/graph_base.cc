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

#include "graph/graph_base.h"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "utils/file_utils.h"

namespace circinus {

std::vector<LabelID> GraphBase::loadUndirectedGraph(const std::string& path) {
  auto infile = openFile(path);

  char line_type;
  // process line: t n_vertices n_edges
  infile >> line_type >> n_vertices_ >> n_edges_;
  vlist_.resize(n_vertices_ + 1);
  vlist_.front() = 0;
  elist_.resize(n_edges_ * 2);
  std::vector<LabelID> labels(n_vertices_);

  // the first n_vertices_ lines should be of type v, with continuous vertex id from 0
  VertexID id;
  LabelID label;
  VertexID degree;
  for (uint32_t i = 0; i < n_vertices_; ++i) {
    CHECK(infile >> line_type >> id >> label >> degree);
    DCHECK_EQ(line_type, 'v');
    DCHECK_EQ(id, i);
    labels[id] = label;
    vlist_[id + 1] = vlist_[id] + degree;
    max_degree_ = std::max(max_degree_, degree);
    vertex_cardinality_by_label_[label] += 1;
  }

  // next n_edges_ lines should be of type e
  std::vector<uint32_t> next_neighbor_offset(n_vertices_, 0);
  VertexID v1, v2;
  for (EdgeID i = 0; i < n_edges_; ++i) {
    CHECK(infile >> line_type) << i + n_vertices_ + 2 << " " << path;
    if (line_type != 'e') {  // there may be a dummy 0 after "e src dst"
      CHECK(infile >> line_type >> v1 >> v2) << i + n_vertices_ + 2 << " " << path;
    } else {
      CHECK(infile >> v1 >> v2) << i + n_vertices_ + 2 << " " << path;
    }
    DCHECK_EQ(line_type, 'e');

    // store an undirected edge as two directed edges
    uint32_t offset = vlist_[v1] + next_neighbor_offset[v1];
    elist_[offset] = v2;
    offset = vlist_[v2] + next_neighbor_offset[v2];
    elist_[offset] = v1;
    ++next_neighbor_offset[v1];
    ++next_neighbor_offset[v2];
  }

  infile.close();

  // check
  for (VertexID i = 0; i < n_vertices_; ++i) {
    DCHECK_EQ(next_neighbor_offset[v1], getVertexOutDegree(v1)) << v1;
  }

  // sort neighbors by id for each vertex
  for (VertexID i = 0; i < n_vertices_; ++i) {
    std::sort(elist_.begin() + vlist_[i], elist_.begin() + vlist_[i + 1]);
  }
  return labels;
}

void GraphBase::loadUndirectedGraphFromBinary(std::istream& input) {
  clear();
  input.read(reinterpret_cast<char*>(&n_vertices_), sizeof(n_vertices_));
  input.read(reinterpret_cast<char*>(&n_edges_), sizeof(n_edges_));
  input.read(reinterpret_cast<char*>(&max_degree_), sizeof(max_degree_));
  binaryStreamToVector(input, vlist_);
  CHECK_EQ(vlist_.size(), n_vertices_ + 1);
  binaryStreamToVector(input, elist_);
  CHECK_EQ(elist_.size(), 2 * n_edges_);
  {  // vertex_cardinality_by_label_
    size_t map_size;
    input.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (size_t i = 0; i < map_size; ++i) {
      LabelID l;
      uint32_t count;
      input.read(reinterpret_cast<char*>(&l), sizeof(l));
      input.read(reinterpret_cast<char*>(&count), sizeof(count));
      vertex_cardinality_by_label_[l] = count;
    }
  }
}

void GraphBase::saveAsBinaryInner(std::ostream& output) const {
  output.write(reinterpret_cast<const char*>(&n_vertices_), sizeof(n_vertices_));
  output.write(reinterpret_cast<const char*>(&n_edges_), sizeof(n_edges_));
  output.write(reinterpret_cast<const char*>(&max_degree_), sizeof(max_degree_));
  vectorToBinaryStream(output, vlist_);
  vectorToBinaryStream(output, elist_);
  {  // vertex_cardinality_by_label_
    auto size = vertex_cardinality_by_label_.size();
    output.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (auto& pair : vertex_cardinality_by_label_) {
      output.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
      output.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second));
    }
  }
}

std::pair<double, double> GraphBase::getMemoryUsage() const {
  std::pair<double, double> ret{0, 0};
  ret.first += sizeof(GraphBase);
  // calculate for unordered_map following the memory usage equation from phmap README
  auto pair_size = sizeof(unordered_map<LabelID, uint32_t>::value_type);
  ret.first += vertex_cardinality_by_label_.bucket_count() * (pair_size + 1);
  ret.second = ret.first;  // for primitive members and hashmap, the two are the same
  ret.first += vlist_.capacity() * sizeof(EdgeID);
  ret.first += elist_.capacity() * sizeof(VertexID);
  ret.second += vlist_.size() * sizeof(EdgeID);
  ret.second += elist_.size() * sizeof(VertexID);
  return ret;
}

}  // namespace circinus

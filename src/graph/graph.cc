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

#include "graph/graph.h"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

namespace circinus {

Graph::Graph(const std::string& path) {
  std::ifstream infile(path);
  CHECK(infile.is_open()) << "Cannot open file " << path;

  char line_type;
  // process line: t n_vertices n_edges
  infile >> line_type >> n_vertices_ >> n_edges_;
  vlist_.resize(n_vertices_ + 1);
  vlist_.front() = 0;
  elist_.resize(n_edges_ * 2);
  labels_.resize(n_vertices_);

  // the first n_vertices_ lines should be of type v, with continuous vertex id from 0
  VertexID id;
  LabelID label;
  VertexID degree;
  for (uint32_t i = 0; i < n_vertices_; ++i) {
    CHECK(infile >> line_type >> id >> label >> degree);
    DCHECK_EQ(line_type, 'v');
    DCHECK_EQ(id, i);
    labels_[id] = label;
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

  // build vertex by label index
  for (auto& pair : vertex_cardinality_by_label_) {
    vertex_ids_by_label_[pair.first].reserve(pair.second);
  }
  // sort neighbors by id for each vertex
  for (VertexID i = 0; i < n_vertices_; ++i) {
    std::sort(elist_.begin() + vlist_[i], elist_.begin() + vlist_[i + 1]);
    vertex_ids_by_label_[labels_[i]].push_back(i);
  }
}

void Graph::DumpToFile(const std::string& path) const {
  std::ofstream output(path);
  output << "t " << getNumVertices() << ' ' << getNumEdges() << '\n';
  for (uint32_t i = 0; i < getNumVertices(); ++i) {
    output << "v " << i << ' ' << getVertexLabel(i) << ' ' << getVertexOutDegree(i) << '\n';
  }

  for (uint32_t i = 0; i < getNumVertices(); ++i) {
    auto neighbors = getOutNeighbors(i);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (i < neighbors.first[j]) {
        output << "e " << i << ' ' << neighbors.first[j] << '\n';
      }
    }
  }
  output.close();
}

void Graph::loadCompressed(std::istream& input) {
  clear();
  input.read(reinterpret_cast<char*>(&n_vertices_), sizeof(n_vertices_));
  input.read(reinterpret_cast<char*>(&n_edges_), sizeof(n_edges_));
  input.read(reinterpret_cast<char*>(&max_degree_), sizeof(max_degree_));
  // vlist
  size_t vlist_size;
  input.read(reinterpret_cast<char*>(&vlist_size), sizeof(vlist_size));
  CHECK_EQ(vlist_size, n_vertices_ + 1);
  vlist_.resize(vlist_size);
  for (size_t i = 0; i < vlist_size; ++i) {
    input.read(reinterpret_cast<char*>(&vlist_[i]), sizeof(vlist_[i]));
  }
  // elist
  size_t elist_size;
  input.read(reinterpret_cast<char*>(&elist_size), sizeof(elist_size));
  CHECK_EQ(elist_size, 2 * n_edges_);
  elist_.resize(elist_size);
  for (size_t i = 0; i < elist_size; ++i) {
    input.read(reinterpret_cast<char*>(&elist_[i]), sizeof(elist_[i]));
  }
  // labels
  size_t labels_size;
  input.read(reinterpret_cast<char*>(&labels_size), sizeof(labels_size));
  CHECK_EQ(labels_size, n_vertices_);
  labels_.resize(labels_size);
  for (size_t i = 0; i < labels_size; ++i) {
    input.read(reinterpret_cast<char*>(&labels_[i]), sizeof(labels_[i]));
  }
  // vertex_cardinality_by_label
  size_t map_size;
  input.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
  for (size_t i = 0; i < map_size; ++i) {
    LabelID l;
    uint32_t count;
    input.read(reinterpret_cast<char*>(&l), sizeof(l));
    input.read(reinterpret_cast<char*>(&count), sizeof(count));
    vertex_cardinality_by_label_[l] = count;
  }
  // build vertex by label index
  for (auto& pair : vertex_cardinality_by_label_) {
    vertex_ids_by_label_[pair.first].reserve(pair.second);
  }
  for (VertexID i = 0; i < n_vertices_; ++i) {
    vertex_ids_by_label_[labels_[i]].push_back(i);
  }
}

void Graph::saveAsBinary(const std::string& path) const {
  std::ofstream output(path, std::ios::out | std::ios::binary);
  output.write(reinterpret_cast<const char*>(&n_vertices_), sizeof(n_vertices_));
  output.write(reinterpret_cast<const char*>(&n_edges_), sizeof(n_edges_));
  output.write(reinterpret_cast<const char*>(&max_degree_), sizeof(max_degree_));
  {
    size_t vlist_size = vlist_.size();
    output.write(reinterpret_cast<const char*>(&vlist_size), sizeof(vlist_size));
    for (auto v : vlist_) {
      output.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
  }
  {
    size_t elist_size = elist_.size();
    output.write(reinterpret_cast<const char*>(&elist_size), sizeof(elist_size));
    for (auto v : elist_) {
      output.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
  }
  {
    auto labels_size = labels_.size();
    output.write(reinterpret_cast<const char*>(&labels_size), sizeof(labels_size));
    for (auto l : labels_) {
      output.write(reinterpret_cast<const char*>(&l), sizeof(l));
    }
  }
  {
    auto size = vertex_cardinality_by_label_.size();
    output.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (auto& pair : vertex_cardinality_by_label_) {
      output.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
      output.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second));
    }
  }
  output.close();
}

}  // namespace circinus

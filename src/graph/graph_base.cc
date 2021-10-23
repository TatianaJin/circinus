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
#include <cmath>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "graph/partitioned_graph.h"
#include "utils/file_utils.h"

namespace circinus {

std::vector<LabelID> GraphBase::loadUndirectedGraphFromEdgeList(const std::string& path) {
  auto infile = openFile(path);

  // process line: t n_vertices n_edges
  std::string line;
  for (uint32_t i = 0; i < 4; ++i) {
    std::getline(infile, line);
    if (i == 2) {
      std::stringstream ss(line);
      char pound;
      std::string nodes;
      std::string edges;
      ss >> pound >> nodes >> n_vertices_ >> edges >> n_edges_;
    }
  }
  vlist_.resize(n_vertices_ + 1);
  vlist_.front() = 0;
  elist_.resize(n_edges_ * 2);
  std::vector<LabelID> labels(n_vertices_);

  // the first n_vertices_ lines should be of type v, with continuous vertex id from 0
  std::vector<std::pair<VertexID, VertexID>> edge_list;
  std::vector<VertexID> degrees(n_vertices_, 0);
  std::unordered_map<VertexID, VertexID> vertex_map;
  VertexID id = 0;
  edge_list.reserve(n_edges_);

  std::vector<std::vector<VertexID>> in_adj(n_vertices_);
  std::vector<std::vector<VertexID>> out_adj(n_vertices_);
  std::size_t found = path.find_last_of("/\\");
  std::string path_prefix = path.substr(0, found);
  // for edge property
  std::string edge_property_path = path_prefix + "/grasper/edge_property/part_0";
  auto edge_property_file = openOutputFile(edge_property_path, std::ios::out);
  // for index
  std::string index_path = path_prefix + "/grasper/index/";

  auto edge_label_file = openOutputFile(index_path + "edge_label", std::ios::out);
  edge_label_file << "edge 1";
  edge_label_file.close();
  auto edge_property_index_file = openOutputFile(index_path + "edge_property_index", std::ios::out);
  edge_property_index_file.close();
  auto vtx_property_index_file = openOutputFile(index_path + "vtx_property_index", std::ios::out);
  vtx_property_index_file.close();

  for (uint64_t i = 0; i < n_edges_; ++i) {
    VertexID v1, v2;
    infile >> v1 >> v2;
    if (vertex_map.find(v1) == vertex_map.end()) {
      vertex_map[v1] = id++;
    }
    if (vertex_map.find(v2) == vertex_map.end()) {
      vertex_map[v2] = id++;
    }
    v1 = vertex_map[v1];
    v2 = vertex_map[v2];
    edge_property_file << v1 << " " << v2 << " " << 1 << " []\n";
    edge_list.emplace_back(std::make_pair(v1, v2));
    CHECK_LE(v1, n_vertices_);
    CHECK_LE(v2, n_vertices_);
    in_adj[v2].push_back(v1);
    out_adj[v1].push_back(v2);
    ++degrees[v1];
    ++degrees[v2];
  }
  edge_property_file.close();

  // for grasper's vertices
  std::string vertices_path = path_prefix + "/grasper/vertices/part_0";
  std::ofstream vertices_file;
  vertices_file.open(vertices_path);
  for (uint32_t i = 0; i < n_vertices_; ++i) {
    vertices_file << i + 1 << " " << in_adj[i].size();
    for (auto& v : in_adj[i]) {
      vertices_file << " " << v + 1;
    }
    vertices_file << " " << out_adj[i].size();
    for (auto& v : out_adj[i]) {
      vertices_file << " " << v + 1;
    }
    vertices_file << "\n";
  }
  vertices_file.close();

  uint64_t labels_num = 100;
  // grapser's vtx label
  auto vtx_label_file = openOutputFile(index_path + "vtx_label", std::ios::out);
  for (uint32_t i = 1; i <= labels_num; ++i) {
    vtx_label_file << i << " " << i << "\n";
  }
  vtx_label_file.close();
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> dist(1, labels_num);
  // for grasper's vertices
  std::string vtx_property_path = path_prefix + "/grasper/vtx_property/part_0";
  std::ofstream vtx_property_file;
  vtx_property_file.open(vtx_property_path);

  // for tve
  std::string tve_path = path.substr(0, path.size() - 9) + ".graph.txt";
  auto tve_file = openOutputFile(tve_path, std::ios::out);
  tve_file << "t " << n_vertices_ << " " << n_edges_ << "\n";

  // for vlabel
  std::string vlabel_path = path.substr(0, path.size() - 9) + ".label.txt";
  auto vlabel_file = openOutputFile(vlabel_path, std::ios::out);

  for (uint32_t i = 0; i < n_vertices_; ++i) {
    labels[i] = dist(mt);
    vtx_property_file << i + 1 << " " << labels[i] + 1 << " []\n";
    tve_file << "v " << i << " " << labels[i] << " " << degrees[i] << "\n";
    vlabel_file << i << " " << labels[i] << "\n";
    vlist_[i + 1] = vlist_[i] + degrees[i];
    max_degree_ = std::max(max_degree_, degrees[i]);
    vertex_cardinality_by_label_[labels[i]] += 1;
  }
  vtx_property_file.close();
  vlabel_file.close();

  // next n_edges_ lines should be of type e
  std::vector<uint32_t> next_neighbor_offset(n_vertices_, 0);
  VertexID v1, v2;
  for (EdgeID i = 0; i < n_edges_; ++i) {
    v1 = edge_list[i].first;
    v2 = edge_list[i].second;
    tve_file << "e " << v1 << " " << v2 << " 0\n";
    // store an undirected edge as two directed edges
    uint32_t offset = vlist_[v1] + next_neighbor_offset[v1];
    elist_[offset] = v2;
    offset = vlist_[v2] + next_neighbor_offset[v2];
    elist_[offset] = v1;
    ++next_neighbor_offset[v1];
    ++next_neighbor_offset[v2];
  }
  tve_file.close();
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

std::vector<LabelID> GraphBase::loadTVEUndirectedGraphToGrasperGraph(const std::string& path) {
  auto infile = openFile(path);

  char line_type;
  // process line: t n_vertices n_edges
  infile >> line_type >> n_vertices_ >> n_edges_;
  vlist_.resize(n_vertices_ + 1);
  vlist_.front() = 0;
  elist_.resize(n_edges_ * 2);
  std::vector<LabelID> labels(n_vertices_);

  std::vector<std::vector<VertexID>> in_adj(n_vertices_);
  std::vector<std::vector<VertexID>> out_adj(n_vertices_);
  std::size_t found = path.find_last_of("/\\");
  std::string path_prefix = path.substr(0, found);
  // for edge property
  std::string edge_property_path = path_prefix + "/grasper/edge_property/part_0";
  auto edge_property_file = openOutputFile(edge_property_path, std::ios::out);
  // for index
  std::string index_path = path_prefix + "/grasper/index/";

  auto edge_label_file = openOutputFile(index_path + "edge_label", std::ios::out);
  edge_label_file << "edge 1";
  edge_label_file.close();
  auto edge_property_index_file = openOutputFile(index_path + "edge_property_index", std::ios::out);
  edge_property_index_file.close();
  auto vtx_property_index_file = openOutputFile(index_path + "vtx_property_index", std::ios::out);
  vtx_property_index_file.close();

  // the first n_vertices_ lines should be of type v, with continuous vertex id from 0
  VertexID id;
  LabelID label;
  VertexID degree;
  LabelID max_label = 0;
  // for grasper's vertices
  std::string vtx_property_path = path_prefix + "/grasper/vtx_property/part_0";
  std::ofstream vtx_property_file;
  vtx_property_file.open(vtx_property_path);
  // for vlabel
  std::string vlabel_path = path.substr(0, path.size() - 9) + ".label.txt";
  auto vlabel_file = openOutputFile(vlabel_path, std::ios::out);
  for (uint32_t i = 0; i < n_vertices_; ++i) {
    CHECK(infile >> line_type >> id >> label >> degree);
    DCHECK_EQ(line_type, 'v');
    DCHECK_EQ(id, i);
    vtx_property_file << id + 1 << " " << label + 1 << " []\n";
    vlabel_file << id << " " << label << "\n";
    labels[id] = label;
    vlist_[id + 1] = vlist_[id] + degree;
    max_degree_ = std::max(max_degree_, degree);
    max_label = std::max(max_label, label);
    vertex_cardinality_by_label_[label] += 1;
  }
  vlabel_file.close();

  // grapser's vtx label
  auto vtx_label_file = openOutputFile(index_path + "vtx_label", std::ios::out);
  for (uint32_t i = 1; i <= max_label + 1; ++i) {
    vtx_label_file << i << " " << i << "\n";
  }
  vtx_label_file.close();

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

    edge_property_file << v1 + 1 << " " << v2 + 1 << " " << 1 << " []\n";
    in_adj[v2].push_back(v1);
    out_adj[v1].push_back(v2);

    // store an undirected edge as two directed edges
    uint32_t offset = vlist_[v1] + next_neighbor_offset[v1];
    elist_[offset] = v2;
    offset = vlist_[v2] + next_neighbor_offset[v2];
    elist_[offset] = v1;
    ++next_neighbor_offset[v1];
    ++next_neighbor_offset[v2];
  }

  edge_property_file.close();
  infile.close();

  // for grasper's vertices
  std::string vertices_path = path_prefix + "/grasper/vertices/part_0";
  std::ofstream vertices_file;
  vertices_file.open(vertices_path);
  for (uint32_t i = 0; i < n_vertices_; ++i) {
    vertices_file << i + 1 << " " << in_adj[i].size();
    for (auto& v : in_adj[i]) {
      vertices_file << " " << v + 1;
    }
    vertices_file << " " << out_adj[i].size();
    for (auto& v : out_adj[i]) {
      vertices_file << " " << v + 1;
    }
    vertices_file << "\n";
  }
  vertices_file.close();

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

std::unique_ptr<GraphBase> GraphBase::loadGraphFromBinary(std::istream& input) {
  bool partitioned_graph = false;
  input.read(reinterpret_cast<char*>(&partitioned_graph), sizeof(bool));
  std::unique_ptr<GraphBase> res;
  if (partitioned_graph) {
    res = std::make_unique<ReorderedPartitionedGraph>();
  } else {
    res = std::make_unique<Graph>();
  }
  res->loadUndirectedGraphFromBinary(input);
  return res;
}

void GraphBase::loadUndirectedGraphFromBinary(std::istream& input) {
  clear();
  input.read(reinterpret_cast<char*>(&n_vertices_), sizeof(n_vertices_));
  input.read(reinterpret_cast<char*>(&n_edges_), sizeof(n_edges_));
  input.read(reinterpret_cast<char*>(&max_degree_), sizeof(max_degree_));
  LOG(INFO) << "Vertices " << n_vertices_ << ", edges " << n_edges_ << ", max degree " << max_degree_;
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
  DLOG(INFO) << "Finish Graph Base Load";
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

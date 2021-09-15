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
#include <unordered_map>
#include <vector>

#include "utils/file_utils.h"

namespace circinus {

void Graph::buildLabelIndex() {
  vertex_ids_by_label_.reserve(vertex_cardinality_by_label_.size());
  // build vertex by label index
  for (auto& pair : vertex_cardinality_by_label_) {
    vertex_ids_by_label_[pair.first].reserve(pair.second);
  }
  for (VertexID i = 0; i < getNumVertices(); ++i) {
    vertex_ids_by_label_[labels_[i]].push_back(i);
  }
}

Graph::Graph(const std::string& path, bool build_label_index) {
  labels_ = loadUndirectedGraph(path);
  if (build_label_index) {
    buildLabelIndex();
  }
}

void Graph::loadUndirectedGraphEdgeList(const std::string& path) { labels_ = loadUndirectedGraphFromEdgeList(path); }

void Graph::dumpToFile(const std::string& path) const {
  auto output = openOutputFile(path);
  output << "t " << getNumVertices() << ' ' << getNumEdges() << '\n';
  for (VertexID i = 0; i < getNumVertices(); ++i) {
    output << "v " << i << ' ' << getVertexLabel(i) << ' ' << getVertexOutDegree(i) << '\n';
  }

  for (VertexID i = 0; i < getNumVertices(); ++i) {
    auto neighbors = getOutNeighbors(i);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (i < neighbors.first[j]) {
        output << "e " << i << ' ' << neighbors.first[j] << '\n';
      }
    }
  }
  output.close();
}

void Graph::loadUndirectedGraphFromBinary(std::istream& input) {
  GraphBase::loadUndirectedGraphFromBinary(input);
  google::FlushLogFiles(google::INFO);
  binaryStreamToVector(input, labels_);
  CHECK_EQ(labels_.size(), n_vertices_);
}

void Graph::saveAsBinaryInner(std::ostream& output) const {
  bool partitioned_graph = false;
  output.write(reinterpret_cast<const char*>(&partitioned_graph), sizeof(bool));
  GraphBase::saveAsBinaryInner(output);
  vectorToBinaryStream(output, labels_);
}

std::pair<double, double> Graph::getMemoryUsage() const {
  auto ret = GraphBase::getMemoryUsage();
  ret.first += sizeof(Graph) - sizeof(GraphBase);
  ret.second += sizeof(Graph) - sizeof(GraphBase);
  ret.first += labels_.capacity() * sizeof(LabelID);
  ret.second += labels_.size() * sizeof(LabelID);
  auto pair_size = sizeof(unordered_map<LabelID, std::vector<VertexID>>::value_type);
  ret.first += vertex_ids_by_label_.bucket_count() * (pair_size + 1);
  ret.second += vertex_ids_by_label_.bucket_count() * (pair_size + 1);
  for (auto& pair : vertex_ids_by_label_) {
    ret.first += pair.second.capacity() * sizeof(VertexID);
  }
  for (auto& pair : vertex_ids_by_label_) {
    ret.second += pair.second.size() * sizeof(VertexID);
  }
  return ret;
}

}  // namespace circinus

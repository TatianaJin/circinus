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
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "./metis.h"
#include "graph/types.h"
#include "utils/file_utils.h"
#include "utils/hashmap.h"

namespace circinus {

class GraphBase {
 protected:
  /* metadata */
  VertexID n_vertices_ = 0;
  EdgeID n_edges_ = 0;
  VertexID max_degree_ = 0;
  unordered_map<LabelID, uint32_t> vertex_cardinality_by_label_;

  /* storage */
  std::vector<EdgeID> vlist_;    // size n_vertices_ + 1, { i: the id of the first edge of vertex i }
  std::vector<VertexID> elist_;  // size n_edges_, { i : the destination vertex id of edge i}

 public:
  static std::unique_ptr<GraphBase> loadGraphFromBinary(std::istream& input);

  GraphBase() {}
  virtual ~GraphBase() {}

  inline VertexID getNumVertices() const { return n_vertices_; }
  inline EdgeID getNumEdges() const { return n_edges_; }
  inline VertexID getGraphMaxDegree() const { return max_degree_; }
  inline LabelID getNumLabels() const { return vertex_cardinality_by_label_.size(); }

  inline uint64_t getVertexCardinalityByLabel(LabelID label) const {
    auto pos = vertex_cardinality_by_label_.find(label);
    if (pos != vertex_cardinality_by_label_.end()) {
      return pos->second;
    }
    return 0;
  }

  inline std::vector<LabelID> getLabels() const {
    std::vector<LabelID> labels;
    labels.reserve(vertex_cardinality_by_label_.size());
    for (auto& pair : vertex_cardinality_by_label_) {
      labels.push_back(pair.first);
    }
    return labels;
  }

  inline VertexID getVertexOutDegree(VertexID id) const { return vlist_[id + 1] - vlist_[id]; }

  /** * @returns a pair { starting neighbor pointer, out degree } */
  inline std::pair<const VertexID*, uint32_t> getOutNeighbors(VertexID id) const {
    DCHECK_LT(id, vlist_.size() - 1);
    return std::make_pair(&elist_[vlist_[id]], vlist_[id + 1] - vlist_[id]);
  }

  inline const EdgeID* getVList() const {  // only for test_metis
    return &vlist_[0];
  }

  inline const VertexID* getEList() const {  // only for test_metis
    return &elist_[0];
  }

  /** Use METIS to compute graph partitions.
   * @param nparts The number of desired partitions.
   * @returns The number of edge cuts among partitions.
   */
  std::pair<std::vector<idx_t>, idx_t> getMetisParts(idx_t nparts) const {
    std::vector<idx_t> parts;
    parts.resize(n_vertices_);
    idx_t* xadj = (idx_t*)getVList();
    idx_t* adjncy = (idx_t*)getEList();
    idx_t nvtxs = getNumVertices();
    idx_t ncon = 1;
    idx_t objval;
    idx_t* part = &parts[0];
    METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);
    return std::make_pair(std::move(parts), objval);
  }

  virtual void clear() {
    vlist_.clear();
    elist_.clear();
    vertex_cardinality_by_label_.clear();
    n_vertices_ = 0;
    n_edges_ = 0;
    max_degree_ = 0;
  }

  /* start of persistence functions */

  /**
   * @param path The file to load graph (undirected).
   *
   * File format
   * t n_vertices n_edges
   * v label_id degree     => n_vertices lines, with continuous vertex id from 0
   * e src_id dst_id 0     => n_edges lines, put the end vertex with smaller id first for each edge
   */
  std::vector<LabelID> loadUndirectedGraph(const std::string& path);
  virtual void dumpToFile(const std::string& path) const = 0;
  inline void saveAsBinary(const std::string& path) const {
    auto output = openOutputFile(path, std::ios::out | std::ios::binary);
    saveAsBinaryInner(output);
    output.close();
  }

  /* end of persistence functions */

  /**
   * @returns The memory usage and the minimum memory usage in bytes.
   */
  virtual std::pair<double, double> getMemoryUsage() const;

 protected:
  virtual void loadUndirectedGraphFromBinary(std::istream& input);
  virtual void saveAsBinaryInner(std::ostream& output) const;

  template <typename T>
  static void vectorToBinaryStream(std::ostream& output, const std::vector<T>& vec) {
    size_t size = vec.size();
    output.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& v : vec) {
      output.write(reinterpret_cast<const char*>(&v), sizeof(T));
    }
  }

  template <typename T>
  static void binaryStreamToVector(std::istream& input, std::vector<T>& vec) {
    size_t size;
    input.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);
    for (size_t i = 0; i < size; ++i) {
      input.read(reinterpret_cast<char*>(&vec[i]), sizeof(T));
    }
  }

  void copyMetadata(const GraphBase& src) {
    n_vertices_ = src.n_vertices_;
    n_edges_ = src.n_edges_;
    max_degree_ = src.max_degree_;
    vertex_cardinality_by_label_ = src.vertex_cardinality_by_label_;
  }
};

}  // namespace circinus

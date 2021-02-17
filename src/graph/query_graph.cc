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

#include "graph/query_graph.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "graph/types.h"

namespace circinus {

QueryGraph::QueryGraph(const std::string& path) {
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
  QueryVertexID id;
  LabelID label;
  QueryVertexID degree;
  for (uint32_t i = 0; i < n_vertices_; ++i) {
    CHECK(infile >> line_type >> id >> label >> degree);
    DCHECK_EQ(line_type, 'v');
    DCHECK_EQ(id, i);
    labels_[id] = label;
    vlist_[id + 1] = vlist_[id] + degree;
    max_degree_ = std::max(max_degree_, degree);
  }

  // next n_edges_ lines should be of type e
  std::vector<uint32_t> next_neighbor_offset(n_vertices_, 0);
  QueryVertexID v1, v2;
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

  // sort neighbors by id for each vertex
  for (QueryVertexID i = 0; i < n_vertices_; ++i) {
    std::sort(elist_.begin() + vlist_[i], elist_.begin() + vlist_[i + 1]);
  }
}
}  // namespace circinus

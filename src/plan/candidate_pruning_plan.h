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

#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "graph/types.h"

namespace circinus {

class CandidatePruningPlan {
 public:
  void newLDFScan(const QueryGraph& q) {
    std::vector<std::pair<LabelID, VertexID>> label_degrees;
    label_degrees.reserve(q.getNumVertices());
    for (QueryVertexID i = 0; i < q.getNumVertices(); ++i) {
      label_degrees.emplace_back(q.getVertexLabel(i), q.getVertexOutDegree(i));
    }
    // FIXME
  }

  void newLDFScan(const QueryGraph& q, const std::vector<QueryVertexID>& vertices_to_scan) {
    std::vector<std::pair<LabelID, VertexID>> label_degrees;
    label_degrees.reserve(vertices_to_scan.size());
    for (auto v : vertices_to_scan) {
      label_degrees.emplace_back(q.getVertexLabel(v), q.getVertexOutDegree(v));
    }
    // FIXME
  }

  void newNLFFilter(const QueryGraph& q) {
    // FIXME
  }

 private:
  LDFScan scan_;                          // phase 1 operator 0
  std::vector<Filter*> phase_1_filters_;  // phase 1 operator 1-n
  ExpandOperator expand_;                 // phase 2 operator 0
  std::vector<Filter*> phase2_filters_;   // phase 2 operator 1-n
  // phase 3 operators
};

}  // namespace circinus

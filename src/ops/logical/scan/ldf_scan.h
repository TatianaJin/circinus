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
#include <numeric>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/scans.h"

namespace circinus {

class LogicalLDFScan {
 private:
  struct ScanWorkload {
    VertexID output_cardinality;
    VertexID scan_cardinality;
  };

  std::vector<QueryVertexID> query_vertices_;  // query vertices to scan
  std::vector<LabelID> labels_;
  std::vector<VertexID> out_degrees_;
  std::vector<VertexID> in_degrees_;

 public:
  LogicalLDFScan() {}

  /**
   * @param q An undirected query graph.
   */
  explicit LogicalLDFScan(const QueryGraph& q) {
    query_vertices_.resize(q.getNumVertices());
    std::iota(query_vertices_.begin(), query_vertices_.end(), 0);
    labels_.reserve(q.getNumVertices());
    out_degrees_.reserve(q.getNumVertices());
    for (QueryVertexID i = 0; i < q.getNumVertices(); ++i) {
      labels_.push_back(q.getVertexLabel(i));
      out_degrees_.push_back(q.getVertexOutDegree(i));
    }
  }

  /**
   * @param q A directed query graph.
   */
  explicit LogicalLDFScan(const DirectedQueryGraph& q) : LogicalLDFScan(static_cast<const QueryGraph&>(q)) {
    in_degrees_.reserve(q.getNumVertices());
    for (QueryVertexID i = 0; i < q.getNumVertices(); ++i) {
      in_degrees_.push_back(q.getVertexInDegree(i));
    }
  }

  /**
   * @param q An undirected query graph.
   * @param query_vertices_to_scan The query vertices to scan for candidates.
   */
  LogicalLDFScan(const QueryGraph& q, const std::vector<QueryVertexID>& query_vertices_to_scan)
      : query_vertices_(query_vertices_to_scan) {
    labels_.reserve(query_vertices_to_scan.size());
    out_degrees_.reserve(query_vertices_to_scan.size());
    for (auto v : query_vertices_to_scan) {
      labels_.push_back(q.getVertexLabel(v));
      out_degrees_.push_back(q.getVertexOutDegree(v));
    }
  }

  /**
   * @param q A directed query graph.
   * @param query_vertices_to_scan The query vertices to scan for candidates.
   */
  LogicalLDFScan(const DirectedQueryGraph& q, const std::vector<QueryVertexID>& query_vertices_to_scan)
      : LogicalLDFScan(static_cast<const QueryGraph&>(q), query_vertices_to_scan) {
    in_degrees_.reserve(query_vertices_to_scan.size());
    for (auto v : query_vertices_to_scan) {
      in_degrees_.push_back(q.getVertexInDegree(v));
    }
  }

  std::vector<std::unique_ptr<Scan>> toPhysicalOperators(const GraphMetadata& metadata, ExecutionConfig& exec) const;

  const std::vector<QueryVertexID>& getQueryVertices() const { return query_vertices_; }

 private:
  static void computeParallelPlan(ScanWorkload workload, ExecutionConfig& exec);
};

}  // namespace circinus

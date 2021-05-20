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

#include "ops/scans.h"

#include <numeric>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"

namespace circinus {

/** Scan the vertices using a label index and filter them by degrees. */
template <bool directed>
class LDFScanBase : public Scan {
 private:
  const LabelID label_ = 0;
  const VertexID min_out_degree_ = 0;
  const VertexID min_in_degree_ = 0;

 public:
  LDFScanBase(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf)
      : Scan(conf), label_(label), min_out_degree_(out_d), min_in_degree_(in_d) {}
  LDFScanBase(LabelID label, VertexID out_d, ExecutionConfig& conf) : LDFScanBase(label, out_d, 0, conf) {}

  void scan(const Graph* g, uint32_t batch_size, ScanContext* ctx) const override {
    auto& input = *g->getVerticesByLabel(label_);
    if (ctx->scan_offset >= ctx->scan_end) return;
    uint32_t count = 0;
    while (batch_size > count && ctx->scan_offset < ctx->scan_end) {
      auto v = input[ctx->scan_offset++];
      if (g->getVertexOutDegree(v) >= min_out_degree_ &&
          (!directed || (static_cast<const DirectedGraph*>(g)->getVertexInDegree(v) >= min_in_degree_))) {
        ++count;
        ctx->candidates.push_back(v);
      }
    }
  }
};
using DirectedLDFScan = LDFScanBase<true>;
using LDFScan = LDFScanBase<false>;

/** Scan the whole graph and filter vertices by label and degrees. */
template <bool directed>
class LDFFullScanBase : public Scan {
 private:
  const LabelID label_ = 0;
  const VertexID min_out_degree_ = 0;
  const VertexID min_in_degree_ = 0;

 public:
  LDFFullScanBase(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf)
      : Scan(conf), label_(label), min_out_degree_(out_d), min_in_degree_(in_d) {}
  LDFFullScanBase(LabelID label, VertexID out_d, ExecutionConfig& conf) : LDFFullScanBase(label, out_d, 0, conf) {}

  void scan(const Graph* g, uint32_t batch_size, ScanContext* ctx) const override {
    if (ctx->scan_offset >= ctx->scan_end) return;
    uint32_t count = 0;
    while (batch_size > count && ctx->scan_offset < ctx->scan_end) {
      auto v = ctx->scan_offset++;
      if (g->getVertexLabel(v) == label_ && g->getVertexOutDegree(v) >= min_out_degree_ &&
          (!directed || (static_cast<const DirectedGraph*>(g)->getVertexInDegree(v) >= min_in_degree_))) {
        ++count;
        ctx->candidates.push_back(g->getVertexGlobalId(v));
      }
    }
  }
};
using DirectedLDFFullScan = LDFFullScanBase<true>;
using LDFFullScan = LDFFullScanBase<false>;

/** Scan the whole graph and filter vertices by degree */
template <bool directed>
class DegreeScanBase : public Scan {
  VertexID min_out_degree_;
  VertexID min_in_degree_;

 public:
  DegreeScanBase(VertexID out_d, VertexID in_d, ExecutionConfig& conf)
      : Scan(conf), min_out_degree_(out_d), min_in_degree_(in_d) {}
  DegreeScanBase(VertexID out_d, ExecutionConfig& conf) : DegreeScanBase(out_d, 0, conf) {}

  void scan(const Graph* g, uint32_t batch_size, ScanContext* ctx) const override {
    if (ctx->scan_offset >= ctx->scan_end) return;
    uint32_t count = 0;
    while (batch_size > count && ctx->scan_offset < ctx->scan_end) {
      auto v = ctx->scan_offset++;
      if (g->getVertexOutDegree(v) >= min_out_degree_ &&
          (!directed || (static_cast<const DirectedGraph*>(g)->getVertexInDegree(v) >= min_in_degree_))) {
        ++count;
        ctx->candidates.push_back(g->getVertexGlobalId(v));
      }
    }
  }
};
using DirectedDegreeScan = DegreeScanBase<true>;
using DegreeScan = DegreeScanBase<false>;

// TODO(tatiana): support versions that apply to graphs sorted by label/degree

/* factory functions */

Scan* Scan::newLDFScan(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf,
                       uint32_t label_pruning_method) {
  if (in_d == 0) {
    if (label_pruning_method == 1) {
      return new LDFScan(label, out_d, conf);
    }
    // TODO(tatiana); handle when label_pruning_method == 2
    return new LDFFullScan(label, out_d, conf);
  }
  if (label_pruning_method == 1) {
    return new DirectedLDFScan(label, out_d, in_d, conf);
  }
  // TODO(tatiana); handle when label_pruning_method == 2
  return new DirectedLDFFullScan(label, out_d, in_d, conf);
}

Scan* Scan::newDegreeScan(VertexID out_d, VertexID in_d, ExecutionConfig& conf) {
  if (in_d == 0) {
    return new DegreeScan(out_d, conf);
  }
  return new DirectedDegreeScan(out_d, in_d, conf);
}

}  // namespace circinus

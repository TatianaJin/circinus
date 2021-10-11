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
#include <string>
#include <utility>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/graph_partition.h"
#include "graph/query_graph.h"
#include "graph/types.h"

namespace circinus {

/** Scan the vertices using a label index and filter them by degrees. */
template <bool directed>
class LDFScanBase : public Scan {
 protected:
  const LabelID label_ = 0;
  const VertexID min_out_degree_ = 0;
  const VertexID min_in_degree_ = 0;

 public:
  LDFScanBase(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf,
              std::string&& name = directed ? "DirectedLDFScan" : "LDFScan")
      : Scan(conf, std::move(name)), label_(label), min_out_degree_(out_d), min_in_degree_(in_d) {}
  LDFScanBase(LabelID label, VertexID out_d, ExecutionConfig& conf) : LDFScanBase(label, out_d, 0, conf) {}

  virtual ~LDFScanBase() {}

  void scan(const Graph* g, ScanContext* ctx) const override {
    auto& input = *g->getVerticesByLabel(label_);
    if (ctx->scan_offset >= ctx->scan_end) return;
    if (ctx->vertex_id == ctx->seed.first) {
      if (input[ctx->scan_offset] > ctx->seed.second || input[ctx->scan_end - 1] < ctx->seed.second) {
        return;
      }
      while (ctx->scan_offset < ctx->scan_end) {
        auto v = input[ctx->scan_offset++];
        if (v != ctx->seed.second) continue;
        if (g->getVertexOutDegree(v) >= min_out_degree_ &&
            (!directed || (static_cast<const DirectedGraph*>(g)->getVertexInDegree(v) >= min_in_degree_)) &&
            validate(*g, v)) {
          ctx->candidates.push_back(v);
        }
      }
    }
    while (ctx->scan_offset < ctx->scan_end) {
      auto v = input[ctx->scan_offset++];
      if (g->getVertexOutDegree(v) >= min_out_degree_ &&
          (!directed || (static_cast<const DirectedGraph*>(g)->getVertexInDegree(v) >= min_in_degree_)) &&
          validate(*g, v)) {
        ctx->candidates.push_back(v);
      }
    }
  }

  void scan(const GraphPartition* g, ScanContext* ctx) const override {
    DCHECK(!directed) << "now graph partition only support undirected graph";
    auto range_start = g->getVertexRangeByLabel(label_).first;
    if (ctx->scan_offset >= ctx->scan_end) return;
    if (ctx->vertex_id == ctx->seed.first) {
      if (range_start + ctx->scan_offset > ctx->seed.second || range_start + (ctx->scan_end - 1) < ctx->seed.second) {
        return;
      }
      while (ctx->scan_offset < ctx->scan_end) {
        auto v = range_start + ctx->scan_offset++;
        if (v != ctx->seed.second) continue;
        if (g->getVertexOutDegree(v) >= min_out_degree_ && validate(*g, v)) {
          ctx->candidates.push_back(v);
        }
        break;
      }
    }
    while (ctx->scan_offset < ctx->scan_end) {
      auto v = range_start + ctx->scan_offset++;
      if (g->getVertexOutDegree(v) >= min_out_degree_ && validate(*g, v)) {
        ctx->candidates.push_back(v);
      }
    }
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << name_ << " (label=" << label_ << ", min_out_degree=" << min_out_degree_
       << (directed ? ", min_in_degree_=" + std::to_string(min_in_degree_) : "") << ", input=" << scan_size_
       << ", shards=" << parallelism_ << ", nfilters=" << filters_.size() << ")";
    return ss.str();
  }
};
using DirectedLDFScan = LDFScanBase<true>;
using LDFScan = LDFScanBase<false>;

/** Scan the whole graph and filter vertices by label and degrees. */
template <bool directed>
class LDFFullScanBase : public LDFScanBase<directed> {
 public:
  LDFFullScanBase(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf)
      : LDFScanBase<directed>(label, out_d, in_d, conf, (directed ? "DirectedLDFFullScan" : "LDFFullScan")) {}
  LDFFullScanBase(LabelID label, VertexID out_d, ExecutionConfig& conf) : LDFFullScanBase(label, out_d, 0, conf) {}

  void scan(const Graph* g, ScanContext* ctx) const override {
    if (ctx->scan_offset >= ctx->scan_end) return;
    while (ctx->scan_offset < ctx->scan_end) {
      auto v = g->getVertexGlobalId(ctx->scan_offset++);
      LOG(INFO) << v;
      if (v != ctx->seed.second) continue;
      if (g->getVertexLabel(v) == this->label_ && g->getVertexOutDegree(v) >= this->min_out_degree_ &&
          (!directed || (static_cast<const DirectedGraph*>(g)->getVertexInDegree(v) >= this->min_in_degree_)) &&
          this->validate(*g, v)) {
        ctx->candidates.push_back(v);
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

  void scan(const Graph* g, ScanContext* ctx) const override {
    if (ctx->scan_offset >= ctx->scan_end) return;
    while (ctx->scan_offset < ctx->scan_end) {
      auto v = g->getVertexGlobalId(ctx->scan_offset++);
      if (g->getVertexOutDegree(v) >= min_out_degree_ &&
          (!directed || (static_cast<const DirectedGraph*>(g)->getVertexInDegree(v) >= min_in_degree_)) &&
          validate(*g, v)) {
        ctx->candidates.push_back(v);
      }
    }
  }

  void scan(const GraphPartition* g, ScanContext* ctx) const override {
    if (ctx->scan_offset >= ctx->scan_end) return;
    auto range_start = g->getPartitionOffset();
    while (ctx->scan_offset < ctx->scan_end) {
      auto v = range_start + ctx->scan_offset++;
      if (g->getVertexOutDegree(v) >= min_out_degree_ && validate(*g, v)) {
        ctx->candidates.push_back(v);
      }
    }
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << (directed ? "Directed" : "") << "DegreeScan (min_out_degree=" << min_out_degree_
       << (directed ? ", min_in_degree_=" + std::to_string(min_in_degree_) : "") << ", input=" << scan_size_
       << ", shards=" << parallelism_ << ")";
    return ss.str();
  }
};
using DirectedDegreeScan = DegreeScanBase<true>;
using DegreeScan = DegreeScanBase<false>;

// TODO(tatiana): support versions that apply to graphs sorted by label/degree

/* factory functions */

std::unique_ptr<Scan> Scan::newLDFScan(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf,
                                       uint32_t label_pruning_method) {
  if (in_d == 0) {
    if (label_pruning_method == 1) {
      return std::make_unique<LDFScan>(label, out_d, conf);
    }
    // TODO(tatiana): handle when label_pruning_method == 2
    return std::make_unique<LDFFullScan>(label, out_d, conf);
  }
  if (label_pruning_method == 1) {
    return std::make_unique<DirectedLDFScan>(label, out_d, in_d, conf);
  }
  // TODO(tatiana): handle when label_pruning_method == 2
  return std::make_unique<DirectedLDFFullScan>(label, out_d, in_d, conf);
}

std::unique_ptr<Scan> Scan::newDegreeScan(VertexID out_d, VertexID in_d, ExecutionConfig& conf) {
  if (in_d == 0) {
    return std::make_unique<DegreeScan>(out_d, conf);
  }
  return std::make_unique<DirectedDegreeScan>(out_d, in_d, conf);
}

}  // namespace circinus

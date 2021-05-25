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
#include <string>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/filters.h"
#include "ops/operator.h"

namespace circinus {

struct ScanContext {
  std::vector<VertexID> candidates;  // output
  VertexID scan_offset = -1;
  VertexID scan_end = 0;

  ScanContext(VertexID offset, VertexID end) : scan_offset(offset), scan_end(end) {}
};

/**
 * The Scan operator should be functor-like so as to be shared by multiple threads to run parallel tasks.
 */
class Scan : public Operator {
 protected:
  VertexID scan_size_ = 0;
  std::string name_;
  std::vector<std::unique_ptr<LocalFilter>> filters_;

 public:
  /**
   * @param label_pruning_method The method to filter vertices by label: 2 for finding label range, 1 for building label
   * index, 0 for full scan.
   */
  static std::unique_ptr<Scan> newLDFScan(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf,
                                          uint32_t label_pruning_method);
  static std::unique_ptr<Scan> newDegreeScan(VertexID out_d, VertexID in_d, ExecutionConfig& conf);

  explicit Scan(ExecutionConfig& conf, std::string&& name = "Scan")
      : Operator(conf.getParallelism()), scan_size_(conf.getInputSize()), name_(std::move(name)) {}

  inline void addFilter(std::unique_ptr<LocalFilter>&& filter) { filters_.push_back(std::move(filter)); }

  inline uint32_t getParallelism() const { return parallelism_; }
  inline ScanContext initScanContext(uint32_t task_idx) const {
    DCHECK_LT(task_idx, parallelism_);
    auto chunk_size = scan_size_ / parallelism_;
    CHECK_NE(chunk_size, 0) << "scan_size=" << scan_size_ << ", parallelism=" << parallelism_;
    if (task_idx < scan_size_ % parallelism_) {
      return ScanContext((chunk_size + 1) * task_idx, chunk_size + 1);
    }
    return ScanContext(chunk_size * task_idx + (scan_size_ % parallelism_), chunk_size);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << name_ << " (input=" << scan_size_ << ", shards=" << parallelism_ << ")";
    return ss.str();
  }

  virtual void scan(const Graph* g, ScanContext* ctx) const = 0;

 protected:
  bool validate(const Graph& g, VertexID v) const {
    for (auto& filter : filters_) {
      if (filter->prune(g, v)) {
        return false;
      }
    }
    return true;
  }
};

}  // namespace circinus

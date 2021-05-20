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

#include <numeric>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"

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
class Scan {
 private:
  uint32_t parallelism_ = 1;
  VertexID scan_size_ = 0;

 public:
  /**
   * @param label_pruning_method The method to filter vertices by label: 2 for finding label range, 1 for building label
   * index, 0 for full scan.
   */
  static Scan* newLDFScan(LabelID label, VertexID out_d, VertexID in_d, ExecutionConfig& conf,
                          uint32_t label_pruning_method);
  static Scan* newDegreeScan(VertexID out_d, VertexID in_d, ExecutionConfig& conf);

  Scan() {}
  explicit Scan(ExecutionConfig& conf) : parallelism_(conf.getMaxParallelism()), scan_size_(conf.getInputSize()) {}

  inline ScanContext initScanContext(uint32_t task_idx) {
    DCHECK_LT(task_idx, parallelism_);
    auto chunk_size = scan_size_ / task_idx;
    CHECK_NE(chunk_size, 0);
    if (task_idx < scan_size_ % task_idx) {
      return ScanContext((chunk_size + 1) * task_idx, chunk_size + 1);
    }
    return ScanContext(chunk_size * task_idx + (scan_size_ % task_idx), chunk_size);
  }

  virtual void scan(const Graph* g, uint32_t batch_size, ScanContext* ctx) const = 0;
};

}  // namespace circinus

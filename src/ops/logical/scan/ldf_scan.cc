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

#include "ops/logical/scan/ldf_scan.h"

#include <algorithm>
#include <memory>

#include "graph/graph_metadata.h"
#include "ops/scans.h"

namespace circinus {

static uint32_t CHUNK_SIZE = 1024;

std::vector<std::unique_ptr<Scan>> LogicalLDFScan::toPhysicalOperators(const GraphMetadata& metadata,
                                                                       ExecutionConfig& exec) const {
  auto size = query_vertices_.size();
  std::vector<std::unique_ptr<Scan>> ret;
  ret.reserve(size);
  for (uint32_t i = 0; i < size; ++i) {
    // estimate scan workload
    ScanWorkload workload;
    workload.output_cardinality = metadata.getNumVertices();
    if (metadata.hasLabelFrequency() && labels_[i] != ALL_LABEL) {
      workload.output_cardinality = std::min(workload.output_cardinality, metadata.getLabelFrequency(labels_[i]));
    }
    if (metadata.hasOutDegreeFrequency()) {
      workload.output_cardinality =
          std::min(workload.output_cardinality, metadata.getNumVerticesWithOutDegreeGE(out_degrees_[i]));
    }
    if (metadata.hasInDegreeFrequency()) {
      workload.output_cardinality =
          std::min(workload.output_cardinality, metadata.getNumVerticesWithInDegreeGE(in_degrees_[i]));
    }
    if (workload.output_cardinality == 0) {  // no vertices satisfying the condition, so no need to scan
      DLOG(INFO) << "WARNING: output_cardinality = 0 out of " << metadata.getNumVertices() << " label " << labels_[i];
      ret.push_back(nullptr);
      continue;
    }
    // instantiate physial scan operators
    if (metadata.isLabeled() && labels_[i] != ALL_LABEL) {
      // 2 for finding label range, 1 for building label index, 0 for full scan
      uint32_t label_pruning_method = metadata.isSortedByLabel() ? 2 : metadata.hasLabelIndex();
      workload.scan_cardinality =
          (label_pruning_method == 0) ? metadata.getNumVertices() : metadata.getLabelFrequency(labels_[i]);
      // compute parallel plan
      computeParallelPlan(workload, exec);
      ret.push_back(Scan::newLDFScan(labels_[i], out_degrees_[i], metadata.isDirected() ? in_degrees_[i] : 0, exec,
                                     label_pruning_method));
    } else {
      workload.scan_cardinality = metadata.getNumVertices();
      // compute parallel plan
      computeParallelPlan(workload, exec);
      ret.push_back(Scan::newDegreeScan(out_degrees_[i], metadata.isDirected() ? in_degrees_[i] : 0, exec));
    }
  }
  return ret;
}

void LogicalLDFScan::computeParallelPlan(ScanWorkload workload, ExecutionConfig& exec) {
  exec.setInputSize(workload.scan_cardinality);
  exec.setParallelism((workload.output_cardinality + CHUNK_SIZE - 1) / CHUNK_SIZE);
  // TODO(tatiana): control scan workload per task? Also, consider I/O later.
}

}  // namespace circinus

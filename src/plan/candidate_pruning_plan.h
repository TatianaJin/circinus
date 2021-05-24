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
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/logical/filter/filter.h"
#include "ops/logical/filter/nlf_filter.h"
#include "ops/logical/scan/ldf_scan.h"

namespace circinus {

class CandidatePruningPlan {
 public:
  CandidatePruningPlan() {}

  inline void newLDFScan(const QueryGraph& q) { scan_ = LogicalLDFScan(q); }

  inline void newLDFScan(const QueryGraph& q, const std::vector<QueryVertexID>& vertices_to_scan) {
    scan_ = LogicalLDFScan(q, vertices_to_scan);
  }

  void newNLFFilter(const QueryGraph& q) {
    if (phase_ == 1) {
      phase_1_filters_.push_back(std::make_unique<LogicalNLFFilter>(q, scan_.getQueryVertices()));
    } else {
      // TODO(tatiana)
      // phase_2_filters_.push_back(new LogicalNLFFilter(q, expand_.getQueryVertices()));
    }
  }

  /**
   * The returned scan operators need to be deleted.
   */
  std::vector<std::unique_ptr<Scan>> getScanOperators(GraphMetadata& metadata, ExecutionConfig& exec_conf) {
    auto ret = scan_.toPhysicalOperators(metadata, exec_conf);
    for (auto& logical_filter : phase_1_filters_) {
      auto filters = logical_filter->toPhysicalOperators(metadata, exec_conf);
      DCHECK_EQ(filters.size(), ret.size());
      for (uint32_t i = 0; i < ret.size(); ++i) {
        ret[i]->addFilter(std::move(filters[i]));
      }
    }
    return ret;
  }

  const auto& getPhase1Filters() const { return phase_1_filters_; }
  const auto& getPhase2Filters() const { return phase_2_filters_; }
  // const auto& getPhase3Filters() const { return phase_3_filters_; }

  uint32_t getPhase() const { return phase_; }
  uint32_t completePhase() { return ++phase_; }
  bool isFinished() const { return finished_; }
  void setFinished() { finished_ = true; }

 private:
  uint32_t phase_ = 1;
  bool finished_ = false;
  LogicalLDFScan scan_;                                               // phase 1 scan operator
  std::vector<std::unique_ptr<LogicalLocalFilter>> phase_1_filters_;  // phase 1 operator 1-n
  // LogicalExpandOperator expand_;                             // phase 2 operator 0
  std::vector<std::unique_ptr<LogicalLocalFilter>> phase_2_filters_;  // phase 2 operator 1-n
  // phase 3 operators
};

}  // namespace circinus

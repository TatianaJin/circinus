#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/logical/filter/cfl_filter.h"
#include "ops/logical/filter/daf_filter.h"
#include "ops/logical/filter/filter.h"
#include "ops/logical/filter/gql_filter.h"
#include "ops/logical/filter/nlf_filter.h"
#include "ops/logical/filter/tso_filter.h"
#include "ops/logical/scan/ldf_scan.h"

namespace circinus {

class CandidatePruningPlan {
 private:
  uint32_t phase_ = 1;
  bool finished_ = false;
  bool partition_result_ = false;  // whether to partition the candidate sets
  LogicalLDFScan scan_;            // phase 1 scan operator
  // LogicalExpandOperator expand_;                             // phase 2 expand operator
  std::vector<std::unique_ptr<LogicalLocalFilter>> local_filters_;            // for phase 1 and 2
  std::vector<std::unique_ptr<LogicalNeighborhoodFilter>> neighbor_filters_;  // phase 3 operators

 public:
  CandidatePruningPlan() {}

  inline void newLDFScan(const QueryGraph& q) { scan_ = LogicalLDFScan(q); }

  inline void newLDFScan(const QueryGraph& q, const std::vector<QueryVertexID>& vertices_to_scan) {
    scan_ = LogicalLDFScan(q, vertices_to_scan);
  }

  void newNLFFilter(const QueryGraph& q, bool maximum_neighbor_degree_filter = false) {
    local_filters_.push_back(
        std::make_unique<LogicalNLFFilter>(q, scan_.getQueryVertices(), maximum_neighbor_degree_filter));
  }

  /**
   * The returned scan operators need to be deleted.
   */
  std::vector<std::unique_ptr<Scan>> getScanOperators(const GraphMetadata& metadata, ExecutionConfig& exec_conf) {
    auto ret = scan_.toPhysicalOperators(metadata, exec_conf);
    for (auto& logical_filter : local_filters_) {
      auto filters = logical_filter->toPhysicalOperators(metadata, exec_conf);
      DCHECK_EQ(filters.size(), ret.size()) << scan_.getQueryVertices().size();
      for (uint32_t i = 0; i < ret.size(); ++i) {
        if (ret[i] != nullptr && filters[i] != nullptr) {
          ret[i]->addFilter(std::move(filters[i]));
        }
      }
    }
    return ret;
  }

  uint32_t getPhase() const { return phase_; }
  uint32_t completePhase() {
    local_filters_.clear();
    neighbor_filters_.clear();
    return ++phase_;
  }
  bool isFinished() const { return finished_; }
  void setFinished() { finished_ = true; }
  bool toPartitionResult() const { return partition_result_; }
  void setPartitionResult(bool to_partition) { partition_result_ = to_partition; }

  const auto& getScanQueryVertices() const { return scan_.getQueryVertices(); }

  void newCFLFilter(const QueryGraph* q, const GraphMetadata& metadata, const std::vector<VertexID>& candidate_size,
                    QueryVertexID seed_qv) {
    neighbor_filters_.emplace_back(std::make_unique<LogicalCFLFilter>(metadata, q, candidate_size, seed_qv));
  }

  void newDAFFilter(const QueryGraph* q, const GraphMetadata& metadata, const std::vector<VertexID>& candidate_size) {
    neighbor_filters_.emplace_back(std::make_unique<LogicalDAFFilter>(metadata, q, candidate_size));
  }

  void newTSOFilter(const QueryGraph* q, const GraphMetadata& metadata, const std::vector<VertexID>& candidate_size) {
    neighbor_filters_.emplace_back(std::make_unique<LogicalTSOFilter>(metadata, q, candidate_size));
  }

  void newGQLFilter(const QueryGraph* q) { neighbor_filters_.push_back(std::make_unique<LogicalGQLFilter>(q)); }

  // std::vector<std::unique_ptr<Extender>> getExtender(const GraphMetadata& metadata, ExecutionConfig& exec_conf) {
  //   std::vector<std::unique_ptr<Extender>> ret;
  //   for (auto& filter : neighbor_filters_) {
  //     auto extender = filter->toPhysicalExtenders(metadata, exec_conf);
  //   }
  // }
  /**
   * The returned filters should be executed sequentially, while each filter is executed in parallel.
   */
  std::vector<std::unique_ptr<NeighborhoodFilter>> getFilterOperators(const GraphMetadata& metadata,
                                                                      ExecutionConfig& exec_conf) {
    std::vector<std::unique_ptr<NeighborhoodFilter>> ret;
    for (auto& filter : neighbor_filters_) {
      auto filters = filter->toPhysicalOperators(metadata, exec_conf);
      for (auto& physical_filter : filters) {
        ret.push_back(std::move(physical_filter));
      }
    }
    return ret;
  }
};

}  // namespace circinus

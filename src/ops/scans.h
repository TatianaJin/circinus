#pragma once

#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph.h"
#include "graph/graph_partition.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/filters.h"
#include "ops/filters/local_filter.h"
#include "ops/operator.h"

namespace circinus {

struct ScanContext {
  std::vector<VertexID> candidates;  // output
  VertexID scan_offset = -1;
  VertexID scan_end = 0;
  QueryVertexID vertex_id;
  std::pair<QueryVertexID, VertexID> seed = std::make_pair(DUMMY_QUERY_VERTEX, 0);

  ScanContext(VertexID offset, VertexID end, QueryVertexID vertex, std::pair<QueryVertexID, VertexID> seed_qv_dv)
      : scan_offset(offset), scan_end(end), vertex_id(vertex), seed(seed_qv_dv) {}
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
      : Operator(conf.getParallelism()), scan_size_(conf.getInputSize()), name_(std::move(name)) {
    // DLOG(INFO) << "Scan size " << scan_size_ << " parallelism " << parallelism_;
  }

  virtual ~Scan() {}

  inline void addFilter(std::unique_ptr<LocalFilter>&& filter) { filters_.push_back(std::move(filter)); }

  inline uint32_t getParallelism() const { return parallelism_; }
  inline ScanContext initScanContext(QueryVertexID vertex, uint32_t task_idx,
                                     std::pair<QueryVertexID, VertexID> seed) const {
    DCHECK_LT(task_idx, parallelism_);
    auto chunk_size = scan_size_ / parallelism_;
    CHECK_NE(chunk_size, 0) << "scan_size=" << scan_size_ << ", parallelism=" << parallelism_;
    if (task_idx < scan_size_ % parallelism_) {
      auto offset = (chunk_size + 1) * task_idx;
      return ScanContext(offset, offset + chunk_size + 1, vertex, seed);
    }
    auto offset = chunk_size * task_idx + (scan_size_ % parallelism_);
    return ScanContext(offset, offset + chunk_size, vertex, seed);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << name_ << " (input=" << scan_size_ << ", shards=" << parallelism_ << ")";
    return ss.str();
  }

  virtual void scan(const Graph* g, ScanContext* ctx) const = 0;
  virtual void scan(const GraphPartition* g, ScanContext* ctx) const {
    LOG(FATAL) << getTypename(*this) << " on GraphPartition not supported";
  }

 protected:
  template <typename GraphView>
  bool validate(const GraphView& g, VertexID v) const {
    for (auto& filter : filters_) {
      if (filter->prune(g, v)) {
        return false;
      }
    }
    return true;
  }
};

}  // namespace circinus

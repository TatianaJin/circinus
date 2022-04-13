#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/types.h"

#include "utils/utils.h"

namespace circinus {

/** Thread-safe functor-style filter */
class TargetFilter {
  std::vector<uint32_t> less_than_indices_;
  std::vector<uint32_t> greater_than_indices_;

 public:
  static std::unique_ptr<TargetFilter> newTargetFilter(std::vector<uint32_t>&& lt_conditions,
                                                       std::vector<uint32_t>&& gt_conditions) {
    return std::make_unique<TargetFilter>(std::move(lt_conditions), std::move(gt_conditions));
  }

  TargetFilter(std::vector<uint32_t>&& lt_conditions, std::vector<uint32_t>&& gt_conditions)
      : less_than_indices_(std::move(lt_conditions)), greater_than_indices_(std::move(gt_conditions)) {
    // LOG(INFO) << "lt " << toString(less_than_indices_) << " gt " << toString(greater_than_indices_);
  }

  void filter(std::vector<VertexID>* targets, const CompressedSubgraphs& group) const;

  bool filter(VertexID target, const CompressedSubgraphs& group) const;

  template <typename View>
  void filterView(View& targets, const CompressedSubgraphs& group) const {
    if (less_than_indices_.empty() && greater_than_indices_.empty()) return;
    auto start = targets.begin();
    auto end = targets.end();
    if (!less_than_indices_.empty()) {
      VertexID less_than = std::numeric_limits<VertexID>::max();  // min of all constraints
      for (auto index : less_than_indices_) {
        auto val = group.getKeyVal(index);
        less_than = std::min(less_than, val);
      }
      end = lowerBound(start, end, less_than);
    }
    if (start == end) {
      targets = View();
      return;
    }
    if (!greater_than_indices_.empty()) {
      VertexID greater_than = 0;  // max of all constraints
      for (auto index : greater_than_indices_) {
        auto val = group.getKeyVal(index);
        greater_than = std::max(greater_than, val);
      }
      ++greater_than;  // greater to greater or equal
      start = lowerBound(start, end, greater_than);
    }
    if (start == end) {
      targets = View();
    } else {
      targets = View(start, end);
    }
  }
};  // class TargetFilter

}  // namespace circinus

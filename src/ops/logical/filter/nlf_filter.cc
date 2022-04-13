#include "ops/logical/filter/nlf_filter.h"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "exec/execution_config.h"
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/filters/nlf_filter.h"
#include "ops/logical/filter/filter.h"

namespace circinus {

std::vector<std::unique_ptr<LocalFilter>> LogicalNLFFilter::toPhysicalOperators(const GraphMetadata& metadata,
                                                                                ExecutionConfig& exec) {
  auto size = out_neighbor_label_frequency_.size();
  std::vector<std::unique_ptr<LocalFilter>> ret;
  ret.reserve(size);
  if (metadata.isSortedByLabel()) {
    for (uint32_t i = 0; i < size; ++i) {
      if (out_neighbor_label_frequency_[i].empty()) {
        ret.push_back(nullptr);
      } else {
        ret.push_back(std::make_unique<QuickNLFFilter>(out_neighbor_label_frequency_[i]));
      }
    }
  } else if (maximum_neighbor_degree_filter_) {
    for (uint32_t i = 0; i < size; ++i) {
      if (out_neighbor_label_frequency_[i].empty()) {
        ret.push_back(nullptr);
      } else {
        ret.push_back(
            std::make_unique<MNDNLFFilter>(std::move(out_neighbor_label_frequency_[i]), maximum_neighbor_degrees_[i]));
      }
    }
  } else {
    for (uint32_t i = 0; i < size; ++i) {
      if (out_neighbor_label_frequency_[i].empty()) {
        ret.push_back(nullptr);
      } else {
        ret.push_back(std::make_unique<NLFFilter>(std::move(out_neighbor_label_frequency_[i])));
      }
    }
  }
  if (metadata.isDirected()) {
  }
  return ret;
}

}  // namespace circinus

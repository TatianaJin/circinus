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

#include "ops/logical/filter/nlf_filter.h"

#include <algorithm>
#include <memory>
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
      ret.push_back(std::make_unique<QuickNLFFilter>(std::move(out_neighbor_label_frequency_[i])));
    }
  } else if (maximum_neighbor_degree_filter_) {
    for (uint32_t i = 0; i < size; ++i) {
      ret.push_back(
          std::make_unique<MNDNLFFilter>(std::move(out_neighbor_label_frequency_[i]), maximum_neighbor_degrees_[i]));
    }
  } else {
    for (uint32_t i = 0; i < size; ++i) {
      ret.push_back(std::make_unique<NLFFilter>(std::move(out_neighbor_label_frequency_[i])));
    }
  }
  if (metadata.isDirected()) {
    // TODO(tatiana): directed version
  }
  return ret;
}

}  // namespace circinus

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

#include "ops/filters/nlf_filter.h"

#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"

namespace circinus {

NLFFilter::NLFFilter(const QueryGraph* query_graph, QueryVertexID query_vid)
    : query_graph_(query_graph), query_vid_(query_vid) {
  // build neighbor label frequency map for query_vid
  auto neighbors = query_graph->getOutNeighbors(query_vid);
  for (uint32_t i = 0; i < neighbors.second; ++i) {
    neighbor_label_frequency_[query_graph->getVertexLabel(neighbors.first[i])] += 1;
  }
}

uint32_t NLFFilter::Filter(const Graph& data_graph, const std::vector<VertexID>& candidates,
                           std::vector<VertexID>* output) {
  uint32_t count = 0;
  for (auto candidate : candidates) {
    auto neighbor_label_frequency = neighbor_label_frequency_;
    auto neighbors = data_graph.getOutNeighbors(candidate);
    // check if the neighbors of the candidate satisfy the neighbor label frequency of the query vertex
    for (uint32_t i = 0; i < neighbors.second; ++i) {
      auto label = data_graph.getVertexLabel(neighbors.first[i]);
      auto pos = neighbor_label_frequency.find(label);
      if (pos != neighbor_label_frequency.end()) {
        if (--pos->second == 0) {  // `label` is satisfied
          neighbor_label_frequency.erase(pos);
          if (neighbor_label_frequency.empty()) {  // all labels are satisfied
            ++count;
            output->push_back(candidate);
            break;
          }
        }
      }
    }
  }
  return count;
}

}  // namespace circinus

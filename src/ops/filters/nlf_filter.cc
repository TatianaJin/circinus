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

#include <algorithm>
#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"

namespace circinus {

NLFFilter::NLFFilter(const QueryGraph* query_graph, QueryVertexID query_vid) {
  // build neighbor label frequency map for query_vid
  auto neighbors = query_graph->getOutNeighbors(query_vid);
  for (uint32_t i = 0; i < neighbors.second; ++i) {
    neighbor_label_frequency_[query_graph->getVertexLabel(neighbors.first[i])] += 1;
  }
}

bool NLFFilter::prune(const Graph& data_graph, VertexID candidate) const {
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
          return false;
          break;
        }
      }
    }
  }
  return true;
}

/**
 * @param g An undirected graph with vertices sorted by label (vertices with smaller label id has smaller id).
 */
inline const VertexID* upperBound(const Graph& g, const VertexID* first, const VertexID* last, LabelID label) {
  if (std::distance(first, last) > 32) {  // binary search when the list is long
    return std::upper_bound(first, last, label,
                            [&g](LabelID label, VertexID v) { return label < g.getVertexLabel(v); });
  }
  // linear scan
  for (; first < last; ++first) {
    if (g.getVertexLabel(*first) > label) {
      return first;
    }
  }
  DCHECK_EQ(g.getVertexLabel(*(first - 1)), label);
  return first;
}

bool QuickNLFFilter::prune(const Graph& data_graph, VertexID candidate) const {
  auto neighbors = data_graph.getOutNeighbors(candidate);
  uint32_t n_checked_label = 0;
  // check if the neighbors of the candidate satisfy the neighbor label frequency of the query vertex
  for (uint32_t i = 0; i < neighbors.second;) {
    auto label = data_graph.getVertexLabel(neighbors.first[i]);
    auto next_pos = upperBound(data_graph, neighbors.first + i + 1, neighbors.first + neighbors.second, label);
    auto frequency = next_pos - (neighbors.first + i);
    i += frequency;
    auto pos = neighbor_label_frequency_.find(label);
    if (pos != neighbor_label_frequency_.end()) {
      ++n_checked_label;
      if (pos->second > frequency) {  // if label frequency is not satisfied, prune this candidate
        return true;
      } else if (n_checked_label == neighbor_label_frequency_.size()) {  // all labels have been checked
        return false;
      }
    }
  }
  return true;
}

}  // namespace circinus

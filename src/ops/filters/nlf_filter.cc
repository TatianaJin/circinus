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
        }
      }
    }
  }
  return true;
}

bool NLFFilter::prune(const GraphPartition& data_graph, VertexID candidate) const {
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
        }
      }
    }
  }
  return true;
}

/**
 * @param g An undirected graph with vertices sorted by label (vertices with smaller label id has smaller id).
 */
template <typename GraphView>
inline const VertexID* upperBound(const GraphView& g, const VertexID* first, const VertexID* last, LabelID label) {
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

bool QuickNLFFilter::prune(const GraphPartition& g, VertexID v) const {
  auto neighbors = g.getOutNeighbors(v);
  uint32_t n_checked_label = 0;

  // check if the neighbors of the candidate satisfy the neighbor label frequency of the query vertex
  if (neighbors.first[0] >= g.getPartitionOffset() &&
      neighbors.first[neighbors.second - 1] < g.getPartitionEnd()) {  // if the neighbors are all in the same partition
    for (uint32_t i = 0; i < neighbors.second;) {
      auto label = g.getVertexLabel(neighbors.first[i]);
      auto next_pos = upperBound(g, neighbors.first + i + 1, neighbors.first + neighbors.second, label);
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
  } else {
    auto neighbor_label_frequency = neighbor_label_frequency_;
    for (uint32_t i = 0; i < neighbors.second; ++i) {
      auto label = g.getVertexLabel(neighbors.first[i]);
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
  }
  return true;
}

bool MNDNLFFilter::prune(const Graph& data_graph, VertexID candidate) const {
  auto neighbors = data_graph.getOutNeighbors(candidate);
  // check if the neighbors of the candidate satisfy the neighbor label frequency of the query vertex
  bool prune = true;
  for (uint32_t i = 0; i < neighbors.second; ++i) {
    if (data_graph.getVertexOutDegree(neighbors.first[i]) >= maximum_neighbor_degree_) {  // pass the mnd filter
      prune = false;
      break;
    }
  }
  if (prune) {
    return true;
  }
  return NLFFilter::prune(data_graph, candidate);
}

}  // namespace circinus

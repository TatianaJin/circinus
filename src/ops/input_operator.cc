#include "ops/input_operator.h"

#include "graph/candidate_set_view.h"
#include "graph/compressed_subgraphs.h"
#include "graph/partitioned_graph.h"
#include "ops/filters/filter.h"

namespace circinus {

std::vector<CompressedSubgraphs> InputOperator::getInputs(const void* g,
                                                          const std::vector<CandidateSetView>& candidates) const {
  if (inputs_are_keys_) {
    return std::vector<CompressedSubgraphs>(candidates[starting_vertex_].begin(), candidates[starting_vertex_].end());
  }
  return std::vector<CompressedSubgraphs>({CompressedSubgraphs(
      newVertexSet(std::vector<VertexID>(candidates[starting_vertex_].begin(), candidates[starting_vertex_].end())))});
}

std::vector<CompressedSubgraphs> PartitionedInputOperator::getInputs(
    const void* data_graph, const std::vector<CandidateSetView>& candidates) const {
  if (!qv_pivots_->empty()) {
#ifndef NDEBUG
    std::stringstream ss;
    ss << qv_pivots_->front().second;
    for (auto& pair : *qv_pivots_) {
      ss << "->" << pair.first;
    }
    LOG(INFO) << "Pruning start vertex by partitioning query vertex " << ss.str();
#endif
    auto g = reinterpret_cast<const ReorderedPartitionedGraph*>(data_graph);
    std::vector<VertexID> candidate[2];
    uint32_t cur_idx = 0;
    for (uint32_t i = 1; i < qv_pivots_->size(); ++i) {
      CHECK_EQ((*qv_pivots_)[i].second, (*qv_pivots_)[i - 1].first) << "qv pivots is not a path.";
    }

    candidate[cur_idx].assign(candidates[qv_pivots_->front().second].begin(),
                              candidates[qv_pivots_->front().second].end());
    for (auto& qv_pivot_pair : *qv_pivots_) {
      candidate[cur_idx ^ 1].clear();
      for (auto data_vertex : candidates[qv_pivot_pair.first]) {
        const auto& nbrs = g->getOutNeighbors(data_vertex);
        if (NeighborhoodFilter::intersectionNotNull(
                std::make_pair(candidate[cur_idx].data(), candidate[cur_idx].size()), nbrs)) {
          candidate[cur_idx ^ 1].push_back(data_vertex);
        }
      }
      cur_idx ^= 1;
    }
    if (inputs_are_keys_) {
      return std::vector<CompressedSubgraphs>(candidate[cur_idx].begin(), candidate[cur_idx].end());
    } else {
      return std::vector<CompressedSubgraphs>({CompressedSubgraphs(newVertexSet(candidate[cur_idx]))});
    }
  }
  return InputOperator::getInputs(data_graph, candidates);
}
}  // namespace circinus

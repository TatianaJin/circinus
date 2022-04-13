#include "ops/filters/gql_filter.h"

namespace circinus {

void GQLFilter::filter(const GraphBase* data_graph, std::vector<std::vector<VertexID>>* candidates,
                       FilterContext* ctx) const {
  for (uint32_t i = ctx->offset; i < ctx->end; ++i) {
    VertexID& data_vertex = (*candidates)[query_vertex_][i];
    if (data_vertex == INVALID_VERTEX_ID) {
      continue;
    }

    if (!verify(data_graph, data_vertex, candidates)) {
      data_vertex = INVALID_VERTEX_ID;
    }
  }
}

bool GQLFilter::verify(const GraphBase* data_graph, const VertexID data_vertex,
                       std::vector<std::vector<VertexID>>* candidates) const {
  const auto& query_vertex_neighbors = query_graph_->getOutNeighbors(query_vertex_);
  const auto& data_vertex_neighbors = data_graph->getOutNeighbors(data_vertex);

  std::vector<uint32_t> left_to_right_edge;
  uint32_t* left_to_right_offset = new uint32_t[query_vertex_neighbors.second + 1];
  uint32_t edge_count = 0;

  for (uint32_t i = 0; i < query_vertex_neighbors.second; ++i) {
    QueryVertexID query_vertex_neighbor = query_vertex_neighbors.first[i];
    left_to_right_offset[i] = edge_count;
    const auto& valid_candidates = (*candidates)[query_vertex_neighbor];
    for (uint32_t j = 0; j < data_vertex_neighbors.second; ++j) {
      VertexID data_vertex_neighbor = data_vertex_neighbors.first[j];

      // binary search to check if the neighbor of data_vertex is in the candidate set of the neighbor of query_vertex
      uint32_t idx = std::lower_bound(valid_candidates.begin(), valid_candidates.end(), data_vertex_neighbor) -
                     valid_candidates.begin();
      if (idx < valid_candidates.size() && valid_candidates[idx] == data_vertex_neighbor) {
        left_to_right_edge.emplace_back(j);
        edge_count++;
      }
    }
  }

  left_to_right_offset[query_vertex_neighbors.second] = edge_count;
  bool res = semiperfectBipartiteMatching(left_to_right_offset, left_to_right_edge, query_vertex_neighbors.second,
                                          data_vertex_neighbors.second);
  delete[] left_to_right_offset;
  return res;
}

}  // namespace circinus

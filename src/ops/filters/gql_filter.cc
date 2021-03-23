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

#include "ops/filters/gql_filter.h"

namespace circinus {

GQLFilter::GQLFilter(const QueryGraph* query_graph, QueryVertexID query_vid,
                     std::unordered_map<QueryVertexID, std::unordered_set<VertexID>>* valid_candidates)
    : FilterBase(query_graph, nullptr), query_vid_(query_vid), valid_candidates_(valid_candidates) {}

void GQLFilter::preFilter(const Graph& data_graph, std::vector<VertexID>& candidates) {
  for (auto& data_vertex : candidates) {
    if (!verify(data_graph, data_vertex)) {
      data_vertex = INVALID_VERTEX_ID;
      valid_candidates_->at(query_vid_).erase(data_vertex);
    }
  }
}

uint32_t GQLFilter::Filter(const Graph& data_graph, std::vector<VertexID>& candidates, std::vector<VertexID>* output) {
  for (auto& data_vertex : candidates) {
    if (data_vertex != INVALID_VERTEX_ID && verify(data_graph, data_vertex)) {
      output->emplace_back(data_vertex);
    }
  }
  return output->size();
}

bool GQLFilter::verify(const Graph& data_graph, VertexID data_vertex) {
  const auto& query_vertex_neighbors = query_graph_->getOutNeighbors(query_vid_);
  const auto& data_vertex_neighbors = data_graph.getOutNeighbors(data_vertex);

  std::vector<uint32_t> left_to_right_edge;
  uint32_t* left_to_right_offset = new uint32_t[query_vertex_neighbors.second + 1];
  uint32_t edge_count = 0;
  for (uint32_t i = 0; i < query_vertex_neighbors.second; ++i) {
    QueryVertexID query_vertex_neighbor = query_vertex_neighbors.first[i];
    left_to_right_offset[i] = edge_count;
    for (uint32_t j = 0; j < data_vertex_neighbors.second; ++j) {
      VertexID data_vertex_neighbor = data_vertex_neighbors.first[j];

      if (valid_candidates_->at(query_vertex_neighbor).count(data_vertex_neighbor) != 0) {
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

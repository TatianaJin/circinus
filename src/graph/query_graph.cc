#include "graph/query_graph.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph/types.h"
#include "utils/file_utils.h"
#include "utils/hashmap.h"

namespace circinus {

inline void QueryGraph::readHeader(std::ifstream& infile, bool directed) {
  char line_type;
  // process line: t n_vertices n_edges
  infile >> line_type >> n_vertices_ >> n_edges_;
  vlist_.resize(n_vertices_ + 1);
  vlist_.front() = 0;
  elist_.resize(n_edges_ * (2 - directed));  // multiply by two if undirected
  labels_.resize(n_vertices_);
}

void QueryGraph::readVertices(std::ifstream& infile) {
  char line_type;
  // the first n_vertices_ lines should be of type v, with continuous vertex id from 0
  QueryVertexID id;
  LabelID label;
  QueryVertexID degree;
  uint32_t unlabeled_count = 0;
  for (uint32_t i = 0; i < n_vertices_; ++i) {
    CHECK(infile >> line_type >> id >> label >> degree);
    DCHECK_EQ(line_type, 'v');
    DCHECK_EQ(id, i);
    labels_[id] = label;
    unlabeled_count += (label == ALL_LABEL);
    vlist_[id + 1] = vlist_[id] + degree;
    max_degree_ = std::max(max_degree_, degree);
    vertex_cardinality_by_label_[label] += 1;
  }
  labelling_ = (unlabeled_count == 0) ? Labeled : ((unlabeled_count == n_vertices_) ? Unlabeled : Mix);
}

template <bool directed>
void QueryGraph::readEdges(std::ifstream& infile) {
  char line_type;
  // next n_edges_ lines should be of type e
  std::vector<uint32_t> next_neighbor_offset(n_vertices_, 0);
  QueryVertexID v1, v2;
  for (EdgeID i = 0; i < n_edges_; ++i) {
    CHECK(infile >> line_type) << "line " << i + n_vertices_ + 2;
    if (line_type != 'e') {  // there may be a dummy 0 after "e src dst"
      CHECK(infile >> line_type >> v1 >> v2) << "line " << i + n_vertices_ + 2;
    } else {
      CHECK(infile >> v1 >> v2) << "line " << i + n_vertices_ + 2;
    }
    DCHECK_EQ(line_type, 'e');
    uint32_t offset = vlist_[v1] + next_neighbor_offset[v1];
    elist_[offset] = v2;
    ++next_neighbor_offset[v1];
    // for undirected graph, store an undirected edge as two directed edges
    if
      constexpr(!directed) {
        offset = vlist_[v2] + next_neighbor_offset[v2];
        elist_[offset] = v1;
        ++next_neighbor_offset[v2];
      }
  }
}

QueryGraph::QueryGraph(const std::string& path) {
  auto infile = openFile(path);
  readHeader(infile);
  readVertices(infile);
  readEdges(infile);
  infile.close();

  sortEdges();
}

void QueryGraph::getInducedSubgraphInner(QueryGraph* ret, const std::vector<QueryVertexID>& vertices) const {
  ret->n_vertices_ = vertices.size();
  ret->vlist_.resize(ret->n_vertices_ + 1, 0);
  ret->labels_.resize(ret->n_vertices_ + 1);
  unordered_map<QueryVertexID, QueryVertexID> vset;  // original id : new id
  vset.reserve(ret->n_vertices_);
  for (QueryVertexID i = 0; i < ret->n_vertices_; ++i) {
    vset[vertices[i]] = i;
  }
  for (QueryVertexID i = 0; i < ret->n_vertices_; ++i) {  // for each vertex, update its label and out neighbors
    auto original = vertices[i];
    ret->labels_[i] = getVertexLabel(original);
    auto neighbors = getOutNeighbors(original);
    for (uint32_t j = 0; j < neighbors.second; ++j) {
      if (vset.count(neighbors.first[j])) {
        ret->elist_.push_back(vset[neighbors.first[j]]);
      }
    }
    ret->vlist_[i + 1] = ret->elist_.size();
    // sort neighbors by id
    std::sort(ret->elist_.begin() + ret->vlist_[i], ret->elist_.begin() + ret->vlist_[i + 1]);
    ret->max_degree_ = std::max(ret->max_degree_, ret->getVertexOutDegree(i));
  }
  ret->n_edges_ = ret->elist_.size() / 2;
}

DirectedQueryGraph::DirectedQueryGraph(const std::string& path) {
  {
    auto infile = openFile(path);
    readHeader(infile);
    readVertices(infile);
    readEdges<true>(infile);
  }
  sortEdges();
  constructInEdges();
}

void DirectedQueryGraph::constructInEdges() {
  /* construct indices for in-edges of vertices */
  in_vlist_.resize(n_vertices_ + 1);
  in_vlist_.front() = 0;
  in_elist_.resize(n_edges_);

  // first compute in-degrees
  for (auto dst : elist_) {
    ++in_vlist_[dst + 1];
  }
  // then compute offsets
  for (QueryVertexID i = 0; i < n_vertices_; ++i) {
    max_in_degree_ = std::max(max_in_degree_, (QueryVertexID)in_vlist_[i + 1]);
    in_vlist_[i + 1] += in_vlist_[i];
  }
  std::vector<uint32_t> next_neighbor_offset(n_vertices_, 0);
  for (QueryVertexID i = 0; i < n_vertices_; ++i) {
    auto pair = getOutNeighbors(i);
    for (uint32_t nb_i = 0; nb_i < pair.second; ++nb_i) {
      uint32_t offset = in_vlist_[pair.first[nb_i]] + next_neighbor_offset[pair.first[nb_i]];
      in_elist_[offset] = i;
      ++next_neighbor_offset[pair.first[nb_i]];
    }
  }
}

}  // namespace circinus

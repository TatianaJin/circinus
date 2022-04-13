#include "ops/scans/ldf_scan.h"

#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"

namespace circinus {

LDFScan::LDFScan(const QueryGraph* query_graph, QueryVertexID query_vid, const Graph* data_graph)
    : query_graph_(query_graph),
      query_vid_(query_vid),
      label_(query_graph->getVertexLabel(query_vid)),
      degree_(query_graph->getVertexOutDegree(query_vid)),
      data_graph_(data_graph) {
  candidates_ = data_graph_->getVerticesByLabel(label_);
  if (candidates_ != nullptr) {
    candidate_end_ = candidates_->size();
  }
}

/**
 * @param candidates The output, candidates of the query vertex in data graph.
 * @returns Return false if there is no record, else true.
 */
bool LDFScan::ScanAll(std::vector<VertexID>* candidates) {
  for (; candidate_offset_ < candidate_end_; ++candidate_offset_) {
    auto v = (*candidates_)[candidate_offset_];
    if (data_graph_->getVertexOutDegree(v) >= degree_) {
      candidates->push_back(v);
    }
  }
  return true;
}

/**
 * @param candidates The output, candidates of the query vertex in data graph.
 * @returns The number of records added to candidates
 */
uint32_t LDFScan::Scan(std::vector<VertexID>* candidates, uint32_t batch_size) {
  if (candidate_offset_ >= candidate_end_) return 0;
  uint32_t count = batch_size;
  while (batch_size > 0 && candidate_offset_ < candidate_end_) {
    auto v = (*candidates_)[candidate_offset_++];
    if (data_graph_->getVertexOutDegree(v) >= degree_) {
      --batch_size;
      candidates->push_back(v);
    }
  }
  return count - batch_size;
}

}  // namespace circinus

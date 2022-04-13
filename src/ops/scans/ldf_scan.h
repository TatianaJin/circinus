#pragma once

#include <vector>

#include "graph/graph.h"
#include "graph/query_graph.h"

namespace circinus {

class LDFScan {
 private:
  const QueryGraph* query_graph_;  // the query graph
  QueryVertexID query_vid_;        // id of the query vertex
  LabelID label_;                  // label of the query vertex
  uint32_t degree_;                // degree of the query vertex
  const std::vector<VertexID>* candidates_;
  size_t candidate_offset_ = 0;
  size_t candidate_end_ = 0;
  const Graph* data_graph_;

 public:
  LDFScan(const QueryGraph* query_graph, QueryVertexID query_vid, const Graph* data_graph);

  /**
   * @param candidates The output, candidates of the query vertex in data graph.
   * @returns Return false if there is no record, else true.
   */
  bool ScanAll(std::vector<VertexID>* candidates);

  /**
   * @param candidates The output, candidates of the query vertex in data graph.
   * @returns The number of records added to candidates
   */
  uint32_t Scan(std::vector<VertexID>* candidates, uint32_t batch_size);
};

}  // namespace circinus

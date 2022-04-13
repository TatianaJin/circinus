#pragma once

#include <vector>

#include "graph/graph.h"
#include "graph/graph_partition.h"
#include "graph/types.h"
#include "utils/utils.h"

namespace circinus {

class LocalFilter {
 public:
  virtual ~LocalFilter() {}

  /** @returns The number of records that passed the filter and are added to output */
  uint32_t filter(const Graph& data_graph, const std::vector<VertexID>& candidates,
                  std::vector<VertexID>* output) const {
    uint32_t count = 0;
    for (auto candidate : candidates) {
      if (prune(data_graph, candidate)) {
        continue;
      }
      ++count;
      output->push_back(candidate);
    }
    return count;
  }

  /** @returns True if vertex v is not a valid mapping */
  virtual bool prune(const Graph& data_graph, VertexID v) const = 0;
  virtual bool prune(const GraphPartition& g, VertexID v) const {
    LOG(FATAL) << getTypename(*this) << " on GraphPartition is not supported";
  }
};

}  // namespace circinus

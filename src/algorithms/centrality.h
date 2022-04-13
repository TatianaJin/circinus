#pragma once

#include <queue>
#include <utility>
#include <vector>

#include "graph/query_graph.h"

namespace circinus {
class ClosenessCentrality {
 private:
  const QueryGraph* graph_;

  double computeClosenessCentrality(QueryVertexID qid) {
    std::queue<QueryVertexID> que;
    std::vector<bool> visited(graph_->getNumVertices());
    std::vector<uint32_t> dis(graph_->getNumVertices(), 0);
    que.push(qid);
    visited[qid] = true;
    double ret = 0;
    while (!que.empty()) {
      QueryVertexID top = que.front();
      que.pop();
      auto[nbrs, cnt] = graph_->getOutNeighbors(top);
      for (uint32_t i = 0; i < cnt; ++i) {
        QueryVertexID nbr = nbrs[i];
        if (visited[nbr] == true) {
          continue;
        }
        dis[nbr] = dis[top] + 1;
        ret += dis[nbr];
        visited[nbr] = true;
        que.push(nbr);
      }
    }
    return 1.0 / ret;
  }

 public:
  explicit ClosenessCentrality(const QueryGraph* graph) : graph_(graph) {}

  std::vector<double> getClosenessCentrality() {
    std::vector<double> ret(graph_->getNumVertices());
    for (QueryVertexID qid = 0; qid < graph_->getNumVertices(); ++qid) {
      ret[qid] = computeClosenessCentrality(qid);
    }
    return std::move(ret);
  }
};
}  // namespace circinus

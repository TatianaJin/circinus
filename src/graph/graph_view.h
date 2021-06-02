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

#pragma once

#include <utility>
#include <vector>

#include "graph/types.h"

namespace circinus {

template <typename G>
class GraphView {
  /** Each graph corresponds to a query graph edge in `Expand` */
  std::vector<G*> graphs_;  // owned, need to delete upon destruction

 public:
  explicit GraphView(std::vector<G*>&& graphs) : graphs_(std::move(graphs)) {}

  ~GraphView() {
    for (auto g : graphs_) {
      delete g;
    }
    graphs_.clear();
  }

  /**
   * nbr_label Hint of the label of out neighbors. The results do not guarantee that the hint is used.
   */
  inline VertexID getVertexOutDegreeWithHint(VertexID id, LabelID nbr_label, uint32_t graph_idx) const {
    return graphs_[graph_idx]->getVertexOutDegreeWithHint(id, nbr_label);
  }

  /**
   * nbr_label Hint of the label of out neighbors. The results do not guarantee that the hint is used.
   */
  inline std::pair<const VertexID*, uint32_t> getOutNeighborsWithHint(VertexID id, LabelID nbr_label,
                                                                      uint32_t graph_idx) const {
    return graphs_[graph_idx]->getOutNeighborsWithHint(id, nbr_label);
  }
};

}  // namespace circinus

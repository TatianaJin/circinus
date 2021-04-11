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

#include <vector>

#include "graph/types.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/traverse_operator.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandEdgeOperator : public TraverseOperator {
 protected:
  uint32_t parent_index_;  // index of parent query vertex in the compressed subgraphs
  uint32_t target_index_;  // index of target query vertex in the compressed subgraphs
  QueryVertexID parent_id_;
  QueryVertexID target_id_;

 public:
  /**
   * Create a new TraverseOperator, which needs to be deleted manually, according to whether the parent and target
   * vertex are compression key or not.
   *
   * @param indices The map from the id of a query vertex to its index in the CompressedSubgraphs
   */
  static TraverseOperator* newExpandEdgeKeyToKeyOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold, SubgraphFilter* filter);

  static TraverseOperator* newExpandEdgeKeyToSetOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold, SubgraphFilter* filter);

  static TraverseOperator* newExpandEdgeSetToKeyOperator(QueryVertexID parent_vertex, QueryVertexID target_vertex,
                                                         const unordered_map<QueryVertexID, uint32_t>& indices,
                                                         const std::vector<uint32_t>& same_label_key_indices,
                                                         const std::vector<uint32_t>& same_label_set_indices,
                                                         uint64_t set_pruning_threshold, SubgraphFilter* filter);

  ExpandEdgeOperator(uint32_t parent_index, uint32_t target_index, QueryVertexID parent, QueryVertexID target,
                     const std::vector<uint32_t>& same_label_key_indices,
                     const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                     SubgraphFilter* filter)
      : TraverseOperator(same_label_key_indices, same_label_set_indices, set_pruning_threshold, filter),
        parent_index_(parent_index),
        target_index_(target_index),
        parent_id_(parent),
        target_id_(target) {
    std::stringstream ss;
    for (auto& v : same_label_key_indices_) {
      ss << ' ' << v;
    }
    for (auto& v : same_label_set_indices_) {
      ss << ' ' << v;
    }
    LOG(INFO) << toString() << " same-label vertices to prune" << ss.str();
  }

  virtual ~ExpandEdgeOperator() {}

  // TODO(tatiana): statistics for how much saving is made
 protected:
  inline void toStringInner(std::stringstream& ss) const {
    ss << ' ' << parent_id_ << " -> " << target_id_;
    if (candidates_ != nullptr) ss << " (" << candidates_->size() << ")";
  }
};

}  // namespace circinus

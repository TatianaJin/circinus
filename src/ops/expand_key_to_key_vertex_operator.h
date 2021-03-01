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

#include <string>
#include <vector>

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"

namespace circinus {

class ExpandKeyToKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandKeyToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    uint32_t output_num = 0;
    while (input_index_ < current_inputs_->size()) {
      std::vector<VertexID> new_keys;
      const auto& input = (*current_inputs_)[input_index_++];
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_[parents_[i]];
        uint32_t key_vid = input.getKeyVal(key);
        if (i == 0) {
          intersect(*candidates_, current_data_graph_->getOutNeighbors(key_vid), &new_keys);
        } else {
          intersectInplace(new_keys, current_data_graph_->getOutNeighbors(key_vid), &new_keys);
        }
        if (new_keys.size() == 0) {
          break;
        }
      }
      if (new_keys.size() != 0) {
        for (VertexID new_key : new_keys) {
          if (input.isExisting(new_key)) {
            continue;
          }
          outputs->emplace_back(input, new_key);
          ++output_num;
        }
      }
      if (output_num >= batch_size) {
        break;
      }
      new_keys.clear();
    }
    return output_num;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToKeyVertexOperator";
    toStringInner(ss);
    if (candidates_ != nullptr) ss << " (" << candidates_->size() << ")";
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    auto ret = new ExpandKeyToKeyVertexOperator(parents_, target_vertex_, query_vertex_indices_);
    DCHECK(candidates_ != nullptr);
    ret->candidates_ = candidates_;
    return ret;
  }
};

}  // namespace circinus

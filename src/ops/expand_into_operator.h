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
#include "ops/traverse_operator.h"

namespace circinus {

class ExpandIntoOperator : public TraverseOperator {
 public:
  ExpandIntoOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                     const std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices)
      : parents_(parents), target_vertex_(target_vertex), query_vertex_indices_(query_vertex_indices) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    uint32_t output_num = 0;
    while (input_index_ < current_inputs_->size()) {
      auto input = (*current_inputs_)[input_index_];
      auto key_vertex_id = input.getKeyVal(query_vertex_indices_[target_vertex_]);
      auto key_out_neighbors = current_data_graph_->getOutNeighbors(key_vertex_id);
      bool add = true;
      for (QueryVertexID vid : parents_) {
        std::vector<VertexID> new_set;
        uint32_t id = query_vertex_indices_[vid];
        intersect(*(input.getSet(id)), key_out_neighbors, &new_set);
        if (new_set.size() == 0) {
          add = false;
          break;
        }
        input.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
      }

      if (add) {
        outputs->emplace_back(std::move(input));
        output_num++;
      }
      input_index_++;

      if (output_num >= batch_size) {
        break;
      }
    }
    return output_num;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandIntoOperator";
    for (auto parent : parents_) {
      DCHECK_EQ(query_vertex_indices_.count(parent), 1);
      ss << ' ' << parent;
    }
    DCHECK_EQ(query_vertex_indices_.count(target_vertex_), 1);
    ss << " -> " << target_vertex_;
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    auto ret = new ExpandIntoOperator(parents_, target_vertex_, query_vertex_indices_);
    return ret;
  }

 protected:
  std::unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;
  std::vector<QueryVertexID> parents_;
  QueryVertexID target_vertex_;
};

}  // namespace circinus

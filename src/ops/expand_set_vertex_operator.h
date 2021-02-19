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

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gperftools/profiler.h"
#include "gtest/gtest.h"

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"

namespace circinus {

class ExpandSetVertexOperator : public ExpandVertexOperator {
 public:
  ExpandSetVertexOperator(const QueryGraph* g, std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                       const std::vector<int>& cover_table, std::unordered_map<QueryVertexID, uint32_t>& query_vertex_indices) : ExpandVertexOperator(g, parents, target_vertex, cover_table, query_vertex_indices) {
  }

  void input(const std::vector<CompressedSubgraphs>* inputs, const Graph* data_graph) {
    inputs_ = inputs;
    data_graph_ = data_graph;
    inputs_idx_ = 0;
  }  

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    while (inputs_idx_ < inputs_->size()) {
      const auto& input = (*inputs_)[inputs_idx_];
      std::vector<VertexID> new_set;
      for (uint32_t i = 0; i < parents_.size(); ++i) {
        uint32_t key = query_vertex_indices_[parents_[i]];
        uint32_t key_vid = input.getKeyVal(key);
        LOG(INFO) << candidates_->size() << " " << data_graph_->getOutNeighbors(key_vid).second;
        for (uint32_t i = 0; i < data_graph_->getOutNeighbors(key_vid).second; ++i) {
          LOG(INFO) << *(data_graph_->getOutNeighbors(key_vid).first + i);
        }
        if (i == 0) {
          intersect(*candidates_, data_graph_->getOutNeighbors(key_vid), &new_set);
        } else {
          intersectInplace(new_set, data_graph_->getOutNeighbors(key_vid), &new_set);
        }
        if(new_set.size() == 0) {
          LOG(INFO) << " no result ";
          break;
          
        }
      }
      if (new_set.size() != 0) {
        outputs->emplace_back(input, std::make_shared<std::vector<VertexID>>(new_set));
        output_num++;
      }
      inputs_idx_++;
    }
    return output_num;
  }
 
 protected:
  const std::vector<CompressedSubgraphs>* inputs_;
  const Graph* data_graph_;
  uint32_t inputs_idx_;
};

}  // namespace circinus

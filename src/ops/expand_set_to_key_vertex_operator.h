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
#include <utility>
#include <vector>

#include "graph/query_graph.h"
#include "ops/expand_vertex_operator.h"
#include "utils/hashmap.h"

namespace circinus {

class ExpandSetToKeyVertexOperator : public ExpandVertexOperator {
 public:
  ExpandSetToKeyVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<false>(outputs, batch_size);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<true>(outputs, batch_size);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandSetToKeyVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    auto ret = new ExpandSetToKeyVertexOperator(parents_, target_vertex_, query_vertex_indices_);
    ret->candidates_ = candidates_;
    return ret;
  }

 protected:
  uint32_t profile() {
    auto& input = (*current_inputs_)[input_index_];
    LOG(INFO) << "---------------------";
    for (auto par : parents_) {
      uint32_t loop_num = 0;
      const auto& parent_match = input.getSet(query_vertex_indices_[par]);
      for (VertexID vid : *parent_match) {
        const auto& out_neighbors = current_data_graph_->getOutNeighbors(vid);
        for (uint32_t i = 0; i < out_neighbors.second; ++i) {
          loop_num++;
        }
      }
      LOG(INFO) << "set out_neighbors num " << loop_num;
    }
    LOG(INFO) << "candidates num " << candidates_->size();
    return 0;
  }

  bool isInCandidates(VertexID key) {
    auto lb = std::lower_bound(candidates_->begin(), candidates_->end(), key);
    return lb != candidates_->end() && *lb == key;
  }

  std::pair<uint32_t, uint32_t> getMinimumParent() {
    uint32_t parent = 0, size = 0xFFFFFFFF;
    auto& input = (*current_inputs_)[input_index_];
    for (auto par : parents_) {
      auto current_size = input.getSet(query_vertex_indices_[par])->size();
      if (current_size < size) {
        size = current_size;
        parent = par;
      }
    }
    return std::make_pair(size, parent);
  }

  // TODO(tatiana): see if hard limit on output size is needed
  template <bool profile>
  inline uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) {
    uint32_t output_num = 0;
    while (input_index_ < current_inputs_->size()) {
      auto min_parent_set = getMinimumParent();
      if (min_parent_set.first * 20 < candidates_->size()) {
        DLOG(INFO) << "fromSetNeighborStrategy";
        output_num += fromSetNeighborStrategy<profile>(outputs, min_parent_set.second);
      } else {
        DLOG(INFO) << "fromCandidateStrategy";
        output_num += fromCandidateStrategy<profile>(outputs);
      }
      input_index_++;
      if (output_num >= batch_size) {
        break;
      }
    }
    return output_num;
  }

  template <bool profile>
  uint32_t fromCandidateStrategy(std::vector<CompressedSubgraphs>* outputs) {
    auto& input = (*current_inputs_)[input_index_];
    uint32_t output_num = 0;
    auto key_map = input.getKeyMap();
    for (VertexID key_vertex_id : *candidates_) {
      if (key_map.count(key_vertex_id)) {
        continue;
      }
      auto key_out_neighbors = current_data_graph_->getOutNeighbors(key_vertex_id);
      // TODO(by) hash key_out_neighbors
      CompressedSubgraphs new_output(input, key_vertex_id);
      bool add = true;
      for (uint32_t set_vid : parents_) {
        std::vector<VertexID> new_set;
        uint32_t id = query_vertex_indices_[set_vid];
        intersect(*input.getSet(id), key_out_neighbors, &new_set);
        if
          constexpr(profile) {
            updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.second, new_set.size());
          }
        if (new_set.empty()) {
          add = false;
          break;
        }
        new_output.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
      }

      if (add) {
        outputs->emplace_back(std::move(new_output));
        output_num++;
      }
    }
    return output_num;
  }

  template <bool profile>
  uint32_t fromSetNeighborStrategy(std::vector<CompressedSubgraphs>* outputs, QueryVertexID min_parent) {
    unordered_set<VertexID> visited;
    auto& input = (*current_inputs_)[input_index_];
    uint32_t output_num = 0;
    const auto& parent_match = input.getSet(query_vertex_indices_[min_parent]);
    auto key_map = input.getKeyMap();
    for (VertexID vid : *parent_match) {
      const auto& out_neighbors = current_data_graph_->getOutNeighbors(vid);
      for (uint32_t i = 0; i < out_neighbors.second; ++i) {
        VertexID key_vertex_id = out_neighbors.first[i];
        if (key_map.count(key_vertex_id)) {
          continue;
        }
        if (visited.insert(key_vertex_id).second) {
          if (!isInCandidates(key_vertex_id)) {
            continue;
          }
          auto key_out_neighbors = current_data_graph_->getOutNeighbors(key_vertex_id);
          // TODO(by) hash key_out_neighbors
          CompressedSubgraphs new_output(input, key_vertex_id);
          bool add = true;
          for (uint32_t set_vid : parents_) {
            std::vector<VertexID> new_set;
            uint32_t id = query_vertex_indices_[set_vid];
            intersect(*input.getSet(id), key_out_neighbors, &new_set);
            if
              constexpr(profile) {
                updateIntersectInfo(input.getSet(id)->size() + key_out_neighbors.second, new_set.size());
              }
            if (new_set.empty()) {
              add = false;
              break;
            }
            new_output.UpdateSets(id, std::make_shared<std::vector<VertexID>>(std::move(new_set)));
          }

          if (add) {
            outputs->emplace_back(std::move(new_output));
            output_num++;
          }
        }
      }
    }
    return output_num;
  }
};

}  // namespace circinus

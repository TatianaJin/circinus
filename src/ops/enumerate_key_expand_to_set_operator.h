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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

class EnumerateKeyExpandToSetOperator : public ExpandVertexOperator {
  const std::vector<QueryVertexID> keys_to_enumerate_;
  std::vector<uint32_t> existing_key_parent_indices_;
  std::vector<int> cover_table_;
  unordered_map<QueryVertexID, uint32_t> enumerate_key_old_indices_;
  std::vector<std::pair<uint32_t, int>> set_old_to_new_pos_;
  bool need_new_input_ = true;
#ifndef USE_FILTER
  // for set pruning, the indices of the sets with the same label as the target in the output
  unordered_set<uint32_t> set_indices_;
#endif
  std::vector<int> enumerated_key_pruning_indices_;

  /* transient */
  std::vector<uint32_t> enumerate_key_idx_;                     // size = keys_to_enumerate_.size();
  std::vector<std::vector<VertexID>*> enumerate_key_pos_sets_;  // size = keys_to_enumerate_.size();
  std::vector<std::vector<VertexID>> target_sets_;  // now we store and reuse the intermediate intersection results
  unordered_set<VertexID> existing_vertices_;
  uint32_t n_exceptions_ = 0;
  CompressedSubgraphs output_;

  /* for profiling */
  std::vector<unordered_set<std::string>> parent_tuple_sets_;
  std::vector<VertexID> parent_tuple_;

 public:
  EnumerateKeyExpandToSetOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                  const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
                                  const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
                                  const std::vector<QueryVertexID>& keys_to_enumerate,
                                  const std::vector<int>& cover_table,
                                  const std::array<std::vector<uint32_t>, 2>& same_label_indices,
                                  std::vector<int>&& enumerated_key_pruning_indices, SubgraphFilter* filter);

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size) override {
    return expandInner<QueryType::Execute>(outputs, batch_size);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size,
                                 uint32_t query_type) override {
    if (query_type == 1) {
      return expandInner<QueryType::Profile>(outputs, batch_size);
    } else if (query_type == 2) {
      return expandInner<QueryType::ProfileWithMiniIntersection>(outputs, batch_size);
    }
  }

  std::string toString() const override;

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new EnumerateKeyExpandToSetOperator(*this);
  }

 private:
  template <QueryType>
  uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t batch_size);

  template <QueryType>
  bool expandInner();
};

}  // namespace circinus

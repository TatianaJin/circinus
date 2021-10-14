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

#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/expand_vertex_traverse_context.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

template <typename G, bool intersect_candidates>
class ExpandKeyToSetVertexOperator : public ExpandVertexOperator {
 public:
  ExpandKeyToSetVertexOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                               const unordered_map<QueryVertexID, uint32_t>& query_vertex_indices,
                               const std::vector<uint32_t>& same_label_key_indices,
                               const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                               std::unique_ptr<SubgraphFilter>&& sfilter = nullptr)
      : ExpandVertexOperator(parents, target_vertex, query_vertex_indices, same_label_key_indices,
                             same_label_set_indices, set_pruning_threshold, std::move(sfilter)) {}

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    if (ctx->getQueryType() == QueryType::Profile) return expandInner<QueryType::Profile>(batch_size, ctx);
    if (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect)
      return expandInner<QueryType::ProfileCandidateSIEffect>(batch_size, ctx);
    CHECK(ctx->getQueryType() == QueryType::ProfileWithMiniIntersection) << "Unknown query type "
                                                                         << (uint32_t)ctx->getQueryType();
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, ctx);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandKeyToSetVertexOperator";
    toStringInner(ss);
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first, input_key_size.second + 1};
  }

 protected:
  template <QueryType profile>
  inline uint32_t expandInner(uint32_t batch_size, TraverseContext* base_ctx) const {
    auto ctx = (ExpandVertexTraverseContext*)base_ctx;
    uint32_t output_num = 0;
    for (; output_num < batch_size && ctx->hasNextInput(); ctx->nextInput()) {
      const auto& input = ctx->getCurrentInput();
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      std::vector<VertexID> new_set;
      expandFromParents<G, profile, intersect_candidates>(input, ctx->getDataGraph<G>(), ctx, parent_indices_,
                                                          exceptions, &new_set);
      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
          // consider reuse of partial intersection results at each parent
          if (isProfileWithMiniIntersectionMode(profile)) {
            std::vector<VertexID> parent_tuple(parents_.size());
            for (uint32_t j = 0; j < parent_indices_.size(); ++j) {
              uint32_t key_vid = input.getKeyVal(parent_indices_[j]);
              parent_tuple[j] = key_vid;
              if (j != 0 || intersect_candidates) {
                ctx->updateDistinctSICount(j, parent_tuple, j);
              }
            }
          }
        }
      if (!new_set.empty()) {
#ifdef USE_FILTER
        auto output = ctx->newOutput(input, std::move(new_set));
        if (filter(*output)) {
          ctx->popOutput();
          continue;
        }
#else
        auto output = ctx->newOutput(input, std::move(new_set), same_label_set_indices_, set_pruning_threshold_);
        if (output == nullptr) {
          continue;
        }
#endif
        ++output_num;
      }
    }
    return output_num;
  }
};

}  // namespace circinus

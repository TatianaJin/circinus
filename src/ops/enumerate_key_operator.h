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
#include <array>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/intersect.h"
#include "graph/types.h"
#include "ops/traverse_operator.h"
#include "ops/traverse_operator_utils.h"
#include "ops/types.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

namespace circinus {

class EnumerateKeyTraverseContext : public ExpandVertexTraverseContext {
 private:
  std::vector<uint32_t> enumerate_key_idx_;                              // size = keys_to_enumerate_.size();
  std::vector<const SingleRangeVertexSetView*> enumerate_key_pos_sets_;  // size = keys_to_enumerate_.size();
  CompressedSubgraphs output_;

 public:
  bool need_new_input = true;
  uint32_t n_exceptions = 0;
  unordered_set<VertexID> existing_vertices;

  EnumerateKeyTraverseContext(const CandidateSetView* candidates, const void* data_graph,
                              std::vector<CompressedSubgraphs>* outputs, QueryType profile,
                              uint32_t n_keys_to_enumerate, QueryVertexID output_graph_size,
                              QueryVertexID output_key_size, uint32_t n_parent_qvs = 0)
      : ExpandVertexTraverseContext(candidates, data_graph, outputs, profile, n_parent_qvs),
        enumerate_key_idx_(n_keys_to_enumerate, 0),
        enumerate_key_pos_sets_(n_keys_to_enumerate),
        output_(output_key_size, output_graph_size) {}

  std::unique_ptr<TraverseContext> clone() const override {
    return std::make_unique<EnumerateKeyTraverseContext>(*this);
  }

  inline auto& getOutput() { return output_; }
  inline bool hasKeyToEnumerate(uint32_t enumerate_key_depth) const {
    return enumerate_key_idx_[enumerate_key_depth] < enumerate_key_pos_sets_[enumerate_key_depth]->size();
  }
  inline VertexID currentKeyToEnumerate(uint32_t enumerate_key_depth) const {
    return (*enumerate_key_pos_sets_[enumerate_key_depth])[enumerate_key_idx_[enumerate_key_depth]];
  }
  inline bool isExistingVertex(VertexID v) const { return existing_vertices.count(v) == 1; }

  void setup(const std::vector<QueryVertexID>& keys_to_enumerate,
             const unordered_map<QueryVertexID, uint32_t>& enumerate_key_old_indices,
             const std::vector<std::pair<uint32_t, int>>& set_old_to_new_pos) {
    // reset index
    enumerate_key_idx_[0] = 0;
    // get sets to enumerate as compression key
    auto& input = getCurrentInput();
    uint32_t depth = 0;
    for (auto v : keys_to_enumerate) {
      auto pos = enumerate_key_old_indices.at(v);
      enumerate_key_pos_sets_[depth] = &(*input.getSet(pos));
      ++depth;
    }
    // set existing matched vertices in output
    std::copy(input.getKeys().begin(), input.getKeys().end(), output_.getKeys().begin());
    for (auto& pair : set_old_to_new_pos) {
      output_.UpdateSets(pair.second, input.getSet(pair.first));
    }
    need_new_input = false;
  }

  inline void nextKeyToEnumerate(uint32_t enumerate_key_depth) { ++enumerate_key_idx_[enumerate_key_depth]; }
  inline void resetKeyToEnumerate(uint32_t enumerate_key_depth) { enumerate_key_idx_[enumerate_key_depth] = 0; }

  inline void resetKeyToEnumerateWithConstraint(uint32_t enumerate_key_depth,
                                                const std::vector<uint32_t>& smaller_indices) {
    VertexID least_vid = 0;
    for (uint32_t idx : smaller_indices) {  // larger than the max of the smaller values
      least_vid = std::max(least_vid, (*enumerate_key_pos_sets_[idx])[enumerate_key_idx_[idx]]);
    }
    ++least_vid;  // +1 as lower bound finds the first element larger than or equal to least_vid
    enumerate_key_idx_[enumerate_key_depth] =
        circinus::lowerBound(enumerate_key_pos_sets_[enumerate_key_depth]->begin(),
                             enumerate_key_pos_sets_[enumerate_key_depth]->end(), least_vid) -
        enumerate_key_pos_sets_[enumerate_key_depth]->begin();
  }
};

template <typename G, bool intersect_candidates>
class EnumerateKeyOperator : public TraverseOperator {
  std::vector<QueryVertexID> keys_to_enumerate_;  // parent query vertices whose matches are enumerated as key
  std::vector<int> cover_table_;
  unordered_map<QueryVertexID, uint32_t> enumerate_key_old_indices_;  // input indices of parent query vertices in set
  std::vector<std::pair<uint32_t, int>> set_old_to_new_pos_;  // input-output indices of other query vertices in set
  std::vector<int> enumerated_key_pruning_indices_;  // i, the indices of sets to be pruned by the i-th enumerated key
  std::vector<std::vector<uint32_t>> po_constraint_adj_;  // i: the indices of smaller keys of the i-th key to enumerate

  unordered_map<QueryVertexID, uint32_t> query_vertex_indices_;

 public:
  EnumerateKeyOperator(QueryVertexID target_vertex,
                       const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
                       const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
                       const std::vector<QueryVertexID>& keys_to_enumerate, const std::vector<int>& cover_table,
                       const std::array<std::vector<uint32_t>, 2>& same_label_indices,
                       std::vector<int>&& enumerated_key_pruning_indices, std::unique_ptr<SubgraphFilter>&& sfilter);

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    if (ctx->getQueryType() == QueryType::Profile) return expandInner<QueryType::Profile>(batch_size, ctx);
    if (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect)
      return expandInner<QueryType::ProfileCandidateSIEffect>(batch_size, ctx);
    CHECK(ctx->getQueryType() == QueryType::ProfileWithMiniIntersection) << "unknown query type "
                                                                         << (uint32_t)ctx->getQueryType();
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, ctx);
  }

  bool extend_vertex() const override { return false; }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    uint32_t output_key_num = 0;
    for (auto& pair : query_vertex_indices_) {
      if (cover_table_[pair.first] == 1) {
        ++output_key_num;
      }
    }
    auto ret = std::make_unique<EnumerateKeyTraverseContext>(
        candidates, graph, outputs, profile, keys_to_enumerate_.size(), query_vertex_indices_.size(), output_key_num);
    return ret;
  }

  std::vector<std::unique_ptr<BipartiteGraph>> computeBipartiteGraphs(
      const Graph* g, const std::vector<CandidateSetView>& candidate_sets) override {
    std::vector<std::unique_ptr<BipartiteGraph>> ret;
    return ret;
  }

  std::vector<std::unique_ptr<GraphPartitionBase>> computeGraphPartitions(
      const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const override {
    std::vector<std::unique_ptr<GraphPartitionBase>> ret;
    return ret;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "EnumerateKeyOperator";
    ss << " (keys to enumerate";
    for (auto v : keys_to_enumerate_) {
      ss << ' ' << v;
    }
    ss << ')';
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first + keys_to_enumerate_.size(), input_key_size.second};
  }

  void setPartialOrder(const PartialOrderConstraintMap& conditions,
                       const unordered_map<QueryVertexID, uint32_t>& seen_vertices) override {
    // set constraints related to targets
    TraverseOperator::setPartialOrder(conditions, seen_vertices);
    // set constraints related to keys to enumerate
    auto n_keys = keys_to_enumerate_.size();
    if (n_keys > 1 && !conditions.empty()) {
      // TODO(tatiana): move the following deduplication and ordering logic to AutomorphismCheck
      unordered_map<QueryVertexID,
                    std::pair<uint32_t, std::pair<unordered_set<QueryVertexID>, unordered_set<QueryVertexID>>>>
          smaller_larger_keys;  // key: {old index, {smaller keys, larger keys}}
      std::vector<uint32_t> keys_old_to_new_index(n_keys, UINT32_MAX);
      uint32_t ordered_count = 0;
      /* find keys to enumerate with constraints */
      for (uint32_t i = 0; i < n_keys; ++i) {
        QueryVertexID qv = keys_to_enumerate_[i];
        auto pos = conditions.find(qv);
        if (pos != conditions.end()) {
          smaller_larger_keys[qv].first = i;
          for (auto smaller : pos->second.first) {
            if (enumerate_key_old_indices_.count(smaller)) {  // the other end of constraint is also key to enumerate
              smaller_larger_keys[qv].second.first.insert(smaller);
              smaller_larger_keys[smaller].second.second.insert(qv);
            }
          }
        }
      }
      if (ordered_count == n_keys) {
        return;  // no partial order constraint among keys to enumerate
      }
      po_constraint_adj_.resize(n_keys);
      /* order keys from smaller to larger by partial order */
      std::queue<QueryVertexID> queue;
      unordered_map<QueryVertexID, uint32_t> checked_smaller_count;
      for (auto& p : smaller_larger_keys) {  // find keys without constraint of smaller keys
        if (p.second.second.first.empty()) {
          queue.push(p.first);
        }
      }
      while (!queue.empty()) {
        auto key = queue.front();
        queue.pop();
        auto& entry = smaller_larger_keys.at(key);
        auto old_index = entry.first;
        unordered_set<QueryVertexID> covered_smaller_keys;
        // find all smaller than smaller keys and deduplicate constraints
        for (auto smaller : entry.second.first) {
          auto& smaller_than_smaller = smaller_larger_keys.at(smaller).second.first;
          covered_smaller_keys.insert(smaller_than_smaller.begin(), smaller_than_smaller.end());
        }
        for (auto smaller : entry.second.first) {
          if (covered_smaller_keys.insert(smaller).second) {  // not covered by any other smaller key
            auto smaller_index = keys_old_to_new_index[smaller_larger_keys.at(smaller).first];
            DCHECK_NE(smaller_index, UINT32_MAX) << "smaller keys should have been set in order";
            po_constraint_adj_[ordered_count].push_back(smaller_index);
          }
        }
        entry.second.first.swap(covered_smaller_keys);
        keys_old_to_new_index[old_index] = ordered_count++;
        // check all larger keys and update queue
        for (auto next_key : entry.second.second) {
          if (++checked_smaller_count[next_key] == smaller_larger_keys.at(next_key).second.first.size()) {
            queue.push(next_key);
          }
        }
      }
      DCHECK_EQ(ordered_count, smaller_larger_keys.size());
      std::vector<uint32_t> keys_to_enumerate_order(n_keys);
      std::vector<int> enumerated_key_pruning_indices_by_order(n_keys);
      for (uint32_t old_index = 0; old_index < n_keys; ++old_index) {
        auto new_index = keys_old_to_new_index[old_index];
        if (new_index == UINT32_MAX) {
          new_index = ordered_count++;
        }
        keys_to_enumerate_order[new_index] = keys_to_enumerate_[old_index];
        enumerated_key_pruning_indices_by_order[new_index] = enumerated_key_pruning_indices_[old_index];
      }
      keys_to_enumerate_.swap(keys_to_enumerate_order);
      enumerated_key_pruning_indices_.swap(enumerated_key_pruning_indices_by_order);
    }
  }

 private:
  template <QueryType>
  uint32_t expandInner(uint32_t batch_size, TraverseContext* ctx) const;

  template <QueryType>
  bool expandInner(TraverseContext* ctx) const;
};

template <typename G, bool intersect_candidates>
EnumerateKeyOperator<G, intersect_candidates>::EnumerateKeyOperator(
    QueryVertexID target_vertex, const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
    const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
    const std::vector<QueryVertexID>& keys_to_enumerate, const std::vector<int>& cover_table,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices, std::vector<int>&& enumerated_key_pruning_indices,
    std::unique_ptr<SubgraphFilter>&& subgraph_filter)
    : TraverseOperator(target_vertex, same_label_indices[1], same_label_indices[0], ~0u, std::move(subgraph_filter)),
      keys_to_enumerate_(keys_to_enumerate),
      cover_table_(cover_table),
      enumerated_key_pruning_indices_(std::move(enumerated_key_pruning_indices)),
      query_vertex_indices_(output_query_vertex_indices) {
  CHECK(subgraph_filter_ != nullptr);
  CHECK_EQ(enumerated_key_pruning_indices_.size(), keys_to_enumerate_.size());

  unordered_set<QueryVertexID> keys_to_enumerate_set(keys_to_enumerate.begin(), keys_to_enumerate.end());

  // get index mapping of enumerated key vertex indices and set vertex indices
  for (auto& pair : input_query_vertex_indices) {
    if (cover_table_[pair.first] != 1) {
      auto new_pos = query_vertex_indices_.at(pair.first);
      set_old_to_new_pos_.emplace_back(pair.second, new_pos);
    }
  }
  for (auto v : keys_to_enumerate_) {
    DCHECK(input_query_vertex_indices.count(v));
    enumerate_key_old_indices_[v] = input_query_vertex_indices.at(v);
  }
}

template <typename G, bool intersect_candidates>
template <QueryType profile>
uint32_t EnumerateKeyOperator<G, intersect_candidates>::expandInner(uint32_t batch_size,
                                                                    TraverseContext* base_ctx) const {
  auto ctx = dynamic_cast<EnumerateKeyTraverseContext*>(base_ctx);
  DCHECK(ctx != nullptr) << "Expect pointer to EnumerateKeyTraverseContext but got " << getTypename(*ctx);

  uint32_t n_outputs = 0;
  while (ctx->hasNextInput()) {
    if (ctx->need_new_input) {
      // find next input with non-empty candidate target set
      auto& input = ctx->getCurrentInput();
      ctx->existing_vertices = input.getKeyMap();
      ctx->n_exceptions = ctx->existing_vertices.size();
      if
        constexpr(isProfileMode(profile)) { ctx->total_num_input_subgraphs += 1; }
      DCHECK_EQ(set_old_to_new_pos_.size(), ctx->getOutput().getNumSets());
      ctx->setup(keys_to_enumerate_, enumerate_key_old_indices_, set_old_to_new_pos_);
    }
    const auto& input = ctx->getCurrentInput();
    const uint32_t enumerate_key_size = keys_to_enumerate_.size();
    uint32_t enumerate_key_depth = ctx->existing_vertices.size() - ctx->n_exceptions;
    if (enumerate_key_depth != 0) {  // if continuing from a previously processed input
      DCHECK_EQ(enumerate_key_depth, enumerate_key_size - 1)
          << "existing_vertices_.size()=" << ctx->existing_vertices.size()
          << ", input.getNumKeys()=" << input.getNumKeys() << '/' << ctx->n_exceptions;
    }
    while (true) {
      while (ctx->hasKeyToEnumerate(enumerate_key_depth)) {
        // update target set
        auto key_vid = ctx->currentKeyToEnumerate(enumerate_key_depth);
        if (ctx->isExistingVertex(key_vid)) {  // skip current key to ensure isomorphism
          ctx->nextKeyToEnumerate(enumerate_key_depth);
          continue;
        }

        if (enumerate_key_depth == enumerate_key_size - 1) {  // the last key query vertex to enumerate, ready to output
          auto& output = ctx->copyOutput(ctx->getOutput());
          // set the enumerated keys in the output
          bool skip = false;
          for (uint32_t key_i = 0; key_i < enumerate_key_size; ++key_i) {
            auto key = ctx->currentKeyToEnumerate(key_i);
            output.UpdateKey(input.getNumKeys() + key_i, key);
            if (enumerated_key_pruning_indices_[key_i] != -1) {
              const auto& pruning_sets = *subgraph_filter_->getPruningSets(enumerated_key_pruning_indices_[key_i]);
              unordered_set<uint32_t> indices(pruning_sets.begin(), pruning_sets.end());
              auto thres = subgraph_filter_->getSetPruningThreshold(enumerated_key_pruning_indices_[key_i]);
#ifdef USE_FILTER
              bool recursive_prune = false;  // only prune by key
#else
              bool recursive_prune = true;  // recursively prune
#endif
              if (output.pruneExistingSets(key, indices, thres, recursive_prune)) {
                skip = true;
                break;
              }
            }
          }
          ctx->nextKeyToEnumerate(enumerate_key_depth);

          if (skip) {
            ctx->popOutput();
            continue;
          }

          if (++n_outputs == batch_size) {
            return n_outputs;
          }
        } else {  // dfs next key depth
          ctx->existing_vertices.insert(key_vid);
          ++enumerate_key_depth;

          if (!po_constraint_adj_.empty() &&
              !po_constraint_adj_[enumerate_key_depth].empty()) {  // start from the first qualified match
            ctx->resetKeyToEnumerateWithConstraint(enumerate_key_depth, po_constraint_adj_[enumerate_key_depth]);
          } else {
            ctx->resetKeyToEnumerate(enumerate_key_depth);  // start from the first match in the next depth
          }
        }
      }
      if (enumerate_key_depth == 0) {  // current input is completely processed
        ctx->need_new_input = true;
        ctx->nextInput();
        break;
      }
      --enumerate_key_depth;
      ctx->existing_vertices.erase(ctx->currentKeyToEnumerate(enumerate_key_depth));
      ctx->nextKeyToEnumerate(enumerate_key_depth);
    }
  }
  return n_outputs;
}

}  // namespace circinus

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
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/intersect.h"
#include "graph/types.h"
#include "ops/expand_vertex_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

class EnumerateTraverseContext : public TraverseContext {
 private:
  std::vector<uint32_t> enumerate_key_idx_;                     // size = keys_to_enumerate_.size();
  std::vector<std::vector<VertexID>*> enumerate_key_pos_sets_;  // size = keys_to_enumerate_.size();
  CompressedSubgraphs output_;

  /* for profiling */
  std::vector<unordered_set<std::string>> parent_tuple_sets_;
  std::vector<VertexID> parent_tuple_;

 public:
  bool need_new_input = true;
  uint32_t n_exceptions = 0;
  std::vector<std::vector<VertexID>> target_sets;  // now we store and reuse the intermediate intersection results
  unordered_set<VertexID> existing_vertices;

  EnumerateTraverseContext(const std::vector<CompressedSubgraphs>* inputs, const void* data_graph, uint32_t input_index,
                           uint32_t input_end_index, uint32_t n_keys_to_enumerate, QueryVertexID output_graph_size,
                           QueryVertexID output_key_size, uint32_t n_parent_qvs = 0)
      : TraverseContext(inputs, data_graph, input_index, input_end_index),
        enumerate_key_idx_(n_keys_to_enumerate, 0),
        enumerate_key_pos_sets_(n_keys_to_enumerate),
        output_(output_key_size, output_graph_size),
        parent_tuple_sets_(n_parent_qvs),
        parent_tuple_(n_parent_qvs),
        target_sets(n_keys_to_enumerate + 1) {}

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
      enumerate_key_pos_sets_[depth] = input.getSet(pos).get();
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

  template <QueryType profile>
  inline void updateIntersection(uint32_t input_size, uint32_t output_size, uint32_t pidx, VertexID key_vid) {
    if
      constexpr(isProfileMode(profile)) {
        updateIntersectInfo(input_size, output_size);
        if
          constexpr(isProfileWithMiniIntersectionMode(profile)) {
            parent_tuple_[pidx] = key_vid;
            distinct_intersection_count +=
                parent_tuple_sets_[pidx].emplace((char*)parent_tuple_.data(), (pidx + 1) * sizeof(VertexID)).second;
          }
      }
  }
};

template <typename G>
class EnumerateKeyExpandToSetOperator : public ExpandVertexOperator {
  const std::vector<QueryVertexID> keys_to_enumerate_;  // parent query vertices whose matches are enumerated as key
  std::vector<uint32_t> existing_key_parent_indices_;   // input indices of parent query vertices in key
  std::vector<int> cover_table_;
  unordered_map<QueryVertexID, uint32_t> enumerate_key_old_indices_;  // input indices of parent query vertices in set
  std::vector<std::pair<uint32_t, int>> set_old_to_new_pos_;  // input-output indices of other query vertices in set
#ifndef USE_FILTER
  // for set pruning, the indices of the sets with the same label as the target in the output
  unordered_set<uint32_t> set_indices_;
#endif
  std::vector<int> enumerated_key_pruning_indices_;  // i, the indices of sets to be pruned by the i-th enumerated key

 public:
  EnumerateKeyExpandToSetOperator(const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
                                  const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
                                  const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
                                  const std::vector<QueryVertexID>& keys_to_enumerate,
                                  const std::vector<int>& cover_table,
                                  const std::array<std::vector<uint32_t>, 2>& same_label_indices,
                                  std::vector<int>&& enumerated_key_pruning_indices, SubgraphFilter* filter);

  uint32_t expand(uint32_t batch_size, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(batch_size, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t batch_size, TraverseContext* ctx) const override {
    if (ctx->query_type == QueryType::Profile) return expandInner<QueryType::Profile>(batch_size, ctx);
    CHECK(ctx->query_type == QueryType::ProfileWithMiniIntersection) << "unknown query type "
                                                                     << (uint32_t)ctx->query_type;
    return expandInner<QueryType::ProfileWithMiniIntersection>(batch_size, ctx);
  }

  std::unique_ptr<TraverseContext> initTraverseContext(const std::vector<CompressedSubgraphs>* inputs,
                                                       const void* graph, uint32_t input_start, uint32_t input_end,
                                                       QueryType profile) const override {
    auto ret = std::make_unique<EnumerateTraverseContext>(
        inputs, graph, input_start, input_end, keys_to_enumerate_.size(), query_vertex_indices_.size(),
        query_vertex_indices_.size() - 1 - query_vertex_indices_.at(target_vertex_),
        isProfileWithMiniIntersectionMode(profile) ? parents_.size() : 0);
    ret->query_type = profile;
    return ret;
  }

  std::vector<std::unique_ptr<GraphPartitionBase>> computeGraphPartitions(
      const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const override;

  std::string toString() const override {
    std::stringstream ss;
    ss << "EnumerateKeyExpandToSetOperator";
    toStringInner(ss);
    ss << " (keys to enumerate";
    for (auto v : keys_to_enumerate_) {
      ss << ' ' << v;
    }
    ss << ')';
    return ss.str();
  }

 private:
  inline const G* getDataGraph(TraverseContext* ctx) const { return (const G*)(ctx->current_data_graph); }

  template <QueryType>
  uint32_t expandInner(uint32_t batch_size, TraverseContext* ctx) const;

  template <QueryType>
  bool expandInner(TraverseContext* ctx) const;
};

template <typename G>
EnumerateKeyExpandToSetOperator<G>::EnumerateKeyExpandToSetOperator(
    const std::vector<QueryVertexID>& parents, QueryVertexID target_vertex,
    const unordered_map<QueryVertexID, uint32_t>& input_query_vertex_indices,
    const unordered_map<QueryVertexID, uint32_t>& output_query_vertex_indices,
    const std::vector<QueryVertexID>& keys_to_enumerate, const std::vector<int>& cover_table,
    const std::array<std::vector<uint32_t>, 2>& same_label_indices, std::vector<int>&& enumerated_key_pruning_indices,
    SubgraphFilter* subgraph_filter)
    : ExpandVertexOperator(parents, target_vertex, output_query_vertex_indices, same_label_indices[1],
                           same_label_indices[0], ~0u, subgraph_filter),
      keys_to_enumerate_(keys_to_enumerate),
      cover_table_(cover_table),
      enumerated_key_pruning_indices_(std::move(enumerated_key_pruning_indices)) {
  CHECK(subgraph_filter != nullptr);
  CHECK_EQ(enumerated_key_pruning_indices_.size(), keys_to_enumerate_.size());

  unordered_set<QueryVertexID> keys_to_enumerate_set(keys_to_enumerate.begin(), keys_to_enumerate.end());
  // get existing key parent indices
  existing_key_parent_indices_.reserve(parents_.size() - keys_to_enumerate_.size());
  for (auto v : parents) {
    if (keys_to_enumerate_set.count(v) == 0) {
      existing_key_parent_indices_.push_back(query_vertex_indices_[v]);
    }
  }

#ifndef USE_FILTER
  // get target pruning set indices
  if (same_label_set_indices_.size() - keys_to_enumerate_.size() > 0) {
    auto pruning_sets = subgraph_filter_->getPruningSets(0);
    set_pruning_threshold_ = subgraph_filter_->getSetPruningThreshold(0);
    CHECK(pruning_sets != nullptr);
    set_indices_.insert(pruning_sets->begin(), pruning_sets->end());
  }
#endif
  for (auto v : keys_to_enumerate_) {
    auto pos =
        std::find(same_label_set_indices_.begin(), same_label_set_indices_.end(), input_query_vertex_indices.at(v));
    if (pos != same_label_set_indices_.end()) same_label_set_indices_.erase(pos);
  }

  // get index mapping of enumerated key vertex indices and set vertex indices
  uint32_t n_input_keys = 0;
  for (auto& pair : input_query_vertex_indices) {
    if (cover_table_[pair.first] != 1) {
      auto new_pos = query_vertex_indices_.at(pair.first);
      set_old_to_new_pos_.emplace_back(pair.second, new_pos);
    } else {
      n_input_keys += (keys_to_enumerate_set.count(pair.first) == 0);
    }
  }
  for (auto v : keys_to_enumerate_) {
    DCHECK(input_query_vertex_indices.count(v));
    enumerate_key_old_indices_[v] = input_query_vertex_indices.at(v);
    CHECK_EQ(query_vertex_indices_[v], n_input_keys) << v;  // assume contiguous indices after existing keys
    ++n_input_keys;
  }
  // assume target is the last set in output
  CHECK_EQ(query_vertex_indices_[target_vertex_], output_query_vertex_indices.size() - n_input_keys - 1);
}

template <typename G>
std::vector<std::unique_ptr<GraphPartitionBase>> EnumerateKeyExpandToSetOperator<G>::computeGraphPartitions(
    const ReorderedPartitionedGraph* g, const std::vector<CandidateScope>& candidate_scopes) const {
  std::vector<std::unique_ptr<GraphPartitionBase>> ret;
  ret.reserve(parents_.size());
  if (parents_.size() - keys_to_enumerate_.size() > 0) {
    // existing key parents first
    unordered_set<QueryVertexID> keys_to_enumerate_set(keys_to_enumerate_.begin(), keys_to_enumerate_.end());
    for (auto v : parents_) {
      if (keys_to_enumerate_set.count(v) == 0) {
        ret.emplace_back(
            GraphPartitionBase::createGraphPartition(candidate_scopes[v], candidate_scopes[target_vertex_], g));
      }
    }
  }
  // then keys to enumerate
  for (auto parent_vertex : keys_to_enumerate_) {
    ret.emplace_back(
        GraphPartitionBase::createGraphPartition(candidate_scopes[parent_vertex], candidate_scopes[target_vertex_], g));
  }
  return ret;
}

template <typename G>
template <QueryType profile>
uint32_t EnumerateKeyExpandToSetOperator<G>::expandInner(uint32_t batch_size, TraverseContext* base_ctx) const {
  auto ctx = dynamic_cast<EnumerateTraverseContext*>(base_ctx);
  DCHECK(ctx != nullptr) << "Expect pointer to EnumerateTraverseContext but got " << getTypename(*ctx);
  uint32_t n_outputs = 0;
  while (ctx->hasNextInput()) {
    if (ctx->need_new_input) {
      // find next input with non-empty candidate target set
      while (ctx->hasNextInput() && expandInner<profile>(ctx)) {
        if
          constexpr(isProfileMode(profile)) {
            ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
          }
        ctx->nextInput();
      }
      if (!ctx->hasNextInput()) {  // all inputs are consumed
        return n_outputs;
      }
      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
        }
      DCHECK_EQ(set_old_to_new_pos_.size(), ctx->getOutput().getNumSets() - 1);
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
        auto pidx = enumerate_key_depth + existing_key_parent_indices_.size();  // parent index
        DCHECK(ctx->target_sets[enumerate_key_depth + 1].empty()) << enumerate_key_depth;
        auto neighbors = getDataGraph(ctx)->getOutNeighborsWithHint(key_vid, target_label_, pidx);
        if (ctx->target_sets[enumerate_key_depth].empty()) {
          DCHECK_EQ(enumerate_key_depth, 0);
          intersect(*candidates_, neighbors, &ctx->target_sets[enumerate_key_depth + 1], ctx->existing_vertices);
        } else {
          intersect(ctx->target_sets[enumerate_key_depth], neighbors, &ctx->target_sets[enumerate_key_depth + 1],
                    ctx->existing_vertices);
        }
        ctx->updateIntersection<profile>(
            (ctx->target_sets[enumerate_key_depth].empty() ? candidates_->size()
                                                           : ctx->target_sets[enumerate_key_depth].size()) +
                neighbors.size(),
            ctx->target_sets[enumerate_key_depth + 1].size(), pidx, key_vid);
        if (ctx->target_sets[enumerate_key_depth + 1].empty()) {
          ctx->nextKeyToEnumerate(enumerate_key_depth);
          continue;
        }
        if (enumerate_key_depth == enumerate_key_size - 1) {  // the last key query vertex to enumerate, ready to output
          auto& target_set = ctx->target_sets.back();
          auto output = ctx->getOutput();
          // set the enumerated keys in the output
          bool skip = false;
          for (uint32_t key_i = 0; key_i < enumerate_key_size; ++key_i) {
            auto key = ctx->currentKeyToEnumerate(key_i);
            output.UpdateKey(input.getNumKeys() + key_i, key);
            if (enumerated_key_pruning_indices_[key_i] != -1) {
              const auto& pruning_sets = *subgraph_filter_->getPruningSets(enumerated_key_pruning_indices_[key_i]);
              unordered_set<uint32_t> indices(pruning_sets.begin(), pruning_sets.end());
              indices.erase(output.getNumSets() - 1);
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
          if (skip) {
            target_set.clear();
            ctx->nextKeyToEnumerate(enumerate_key_depth);
            continue;
          }
#ifdef USE_FILTER
          // set the target set
          output.UpdateSets(output.getNumSets() - 1, std::make_shared<std::vector<VertexID>>(std::move(target_set)));
          if (filter(output)) {
            ctx->nextKeyToEnumerate(enumerate_key_depth);
            continue;
          }
#else
          // set the target set
          if (target_set.size() == 1) {
            auto indices = set_indices_;
            if (output.pruneExistingSets(target_set.front(), indices, set_pruning_threshold_)) {
              ctx->nextKeyToEnumerate(enumerate_key_depth);
              continue;
            }
          }
          output.UpdateSets(output.getNumSets() - 1, std::make_shared<std::vector<VertexID>>(std::move(target_set)));
#endif

          ctx->outputs->push_back(std::move(output));
          ctx->nextKeyToEnumerate(enumerate_key_depth);
          if (++n_outputs == batch_size) {
            return n_outputs;
          }
        } else {  // dfs next key depth
          ctx->existing_vertices.insert(key_vid);
          ++enumerate_key_depth;
          ctx->resetKeyToEnumerate(enumerate_key_depth);  // start from the first match in the next depth
        }
      }
      ctx->target_sets[enumerate_key_depth].clear();
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

template <typename G>
template <QueryType profile>
bool EnumerateKeyExpandToSetOperator<G>::expandInner(
    TraverseContext* base_ctx) const {  // handles a new input and init the transient states
  auto ctx = (EnumerateTraverseContext*)base_ctx;
  auto& input = ctx->getCurrentInput();
  auto& target_set = ctx->target_sets.front();
  ctx->existing_vertices = input.getKeyMap();
  input.getExceptions(ctx->existing_vertices, {}, same_label_set_indices_);
  ctx->n_exceptions = ctx->existing_vertices.size();

  for (uint32_t i = 0; i < existing_key_parent_indices_.size(); ++i) {
    uint32_t key_vid = input.getKeyVal(existing_key_parent_indices_[i]);
    auto neighbors = getDataGraph(ctx)->getOutNeighborsWithHint(key_vid, target_label_, i);
    auto input_size = target_set.size() + neighbors.size();
    if (i == 0) {
      intersect(*candidates_, neighbors, &target_set, ctx->existing_vertices);
    } else {
      intersectInplace(target_set, neighbors, &target_set);
    }
    ctx->updateIntersection<profile>(input_size, target_set.size(), i, key_vid);
    if (target_set.empty()) {
      return true;
    }
  }

  return candidates_->empty() || (target_set.empty() && existing_key_parent_indices_.size() > 0);
}

}  // namespace circinus

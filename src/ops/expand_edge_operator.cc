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

#include "ops/expand_edge_operator.h"

#include <algorithm>
#include <chrono>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/graph_view.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/traverse_operator.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

#define makeVertexSet(vertex) std::make_shared<std::vector<VertexID>>(std::vector<VertexID>({vertex}))
#define makeShared(set) std::make_shared<std::vector<VertexID>>(std::move(set))

#define CONSTRUCT(name)                                                                                                \
  ExpandEdge##name##Operator(uint32_t parent_index, uint32_t target_index, QueryVertexID parent, QueryVertexID target, \
                             const std::vector<uint32_t>& same_label_key_indices,                                      \
                             const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,      \
                             SubgraphFilter* filter)                                                                   \
      : ExpandEdgeOperator(parent_index, target_index, parent, target, same_label_key_indices, same_label_set_indices, \
                           set_pruning_threshold, filter) {}

template <typename G>
class ExpandEdgeKeyToSetOperator : public ExpandEdgeOperator {
  unordered_set<VertexID> candidate_set_;
  unordered_set<VertexID> parent_set_;  // for profile

 public:
  CONSTRUCT(KeyToSet)

  void setCandidateSets(const std::vector<VertexID>* candidates) override {
    candidates_ = candidates;
    candidate_set_.insert(candidates->begin(), candidates->end());
  }

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    uint32_t n = 0;
    for (; n < cap && input_index_ < current_inputs_->size(); ++input_index_) {
      n += expandInner<QueryType::Execute>(outputs, (*current_inputs_)[input_index_]);
    }
    return n;
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t cap,
                                 uint32_t query_type) override {
    uint32_t old_input_index = input_index_;
    uint32_t n = 0;
    for (; n < cap && input_index_ < current_inputs_->size(); ++input_index_) {
      if (query_type == 1) {
        n += expandInner<QueryType::Profile>(outputs, (*current_inputs_)[input_index_]);
      } else {
        CHECK_EQ(query_type, 2) << "unknown query type " << query_type;
        n += expandInner<QueryType::ProfileWithMiniIntersection>(outputs, (*current_inputs_)[input_index_]);
      }
      total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
    }
    intersection_count_ += input_index_ - old_input_index;
    return n;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandEdgeKeyToSetOperator";
    toStringInner(ss);
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new ExpandEdgeKeyToSetOperator(*this);
  }

 private:
  /** @returns True if one CompressedSubgraphs is generated, else false. */
  template <QueryType profile>
  inline bool expandInner(std::vector<CompressedSubgraphs>* outputs, const CompressedSubgraphs& input) {
    std::vector<VertexID> targets;
    auto parent_match = input.getKeyVal(parent_index_);
    intersect(candidate_set_, ((G*)current_data_graph_)->getOutNeighborsWithHint(parent_match, 0, 0), &targets,
              input.getExceptions(same_label_key_indices_, same_label_set_indices_));
    if
      constexpr(isProfileMode(profile)) {
        if
          constexpr(isProfileWithMiniIntersectionMode(profile)) {
            distinct_intersection_count_ += parent_set_.insert(parent_match).second;
          }
        total_intersection_input_size_ +=
            candidate_set_.size() + ((G*)current_data_graph_)->getVertexOutDegreeWithHint(parent_match, 0, 0);
        total_intersection_output_size_ += targets.size();
      }
    if (targets.empty()) {
      return false;
    }
#ifdef USE_FILTER
    CompressedSubgraphs output(input, std::move(targets));
    if (filter(output)) {  // actively prune existing sets
      return false;
    }
#else
    CompressedSubgraphs output(input, std::move(targets), same_label_set_indices_, set_pruning_threshold_);
    if (output.empty()) {  // actively prune existing sets
      return false;
    }
#endif
    outputs->emplace_back(std::move(output));
    return true;
  }
};

template <typename G, bool intersect_candidates = true>
class ExpandEdgeKeyToKeyOperator : public ExpandEdgeOperator {
  // calculated from current_inputs_[input_index_]
  std::vector<VertexID> current_targets_;
  uint32_t current_target_index_ = 0;
  unordered_set<VertexID> candidate_set_;
  unordered_set<VertexID> parent_set_;  // for profile

 public:
  CONSTRUCT(KeyToKey)

  void setCandidateSets(const std::vector<VertexID>* candidates) override {
    candidates_ = candidates;
    candidate_set_.insert(candidates->begin(), candidates->end());
  }

  template <QueryType profile>
  uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) {
    uint32_t n = 0;
    while (true) {
      // if there are existing targets from the last input, consume first
      if (current_target_index_ < current_targets_.size()) {
        auto& input = (*current_inputs_)[input_index_ - 1];
        while (current_target_index_ < current_targets_.size()) {
#ifdef USE_FILTER
          CompressedSubgraphs output(input, current_targets_[current_target_index_], same_label_set_indices_,
                                     set_pruning_threshold_, false);
#else
          CompressedSubgraphs output(input, current_targets_[current_target_index_], same_label_set_indices_,
                                     set_pruning_threshold_);
#endif
          ++current_target_index_;
          if (output.empty()) continue;
#ifdef USE_FILTER
          if (filter(output)) continue;
#endif
          outputs->emplace_back(std::move(output));
          if (++n == cap) {
            return n;
          }
        }
        if (n == cap) {
          return n;
        }
      }
      // return if all inputs in the current batch are consumed
      if (input_index_ == current_inputs_->size()) {
        return n;
      }
      // consume the next input
      expandInner<profile>((*current_inputs_)[input_index_]);
      if (isProfileMode(profile)) {
        total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
      }
      ++input_index_;
    }
    return n;
  }

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    return expandInner<QueryType::Execute>(outputs, cap);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t cap,
                                 uint32_t query_type) override {
    auto old_input_index = input_index_;
    uint32_t n = 0;
    if (query_type == 1) {
      n = expandInner<QueryType::Profile>(outputs, cap);
    } else {
      CHECK_EQ(query_type, 2) << "unknown query type " << query_type;
      n = expandInner<QueryType::ProfileWithMiniIntersection>(outputs, cap);
    }
    intersection_count_ += (input_index_ - old_input_index) * (intersect_candidates);
    return n;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandEdgeKeyToKeyOperator";
    toStringInner(ss);
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new ExpandEdgeKeyToKeyOperator(*this);
  }

 private:
  template <QueryType profile>
  inline void expandInner(const CompressedSubgraphs& input) {
    current_targets_.clear();
    current_target_index_ = 0;
    auto parent_match = input.getKeyVal(parent_index_);
    if (!intersect_candidates) {
      removeExceptions(((const G*)current_data_graph_)->getOutNeighborsWithHint(parent_match, 0, 0), &current_targets_,
                       input.getExceptions(same_label_key_indices_, same_label_set_indices_));
    } else {
      intersect(candidate_set_, ((const G*)current_data_graph_)->getOutNeighborsWithHint(parent_match, 0, 0),
                &current_targets_, input.getExceptions(same_label_key_indices_, same_label_set_indices_));
      if
        constexpr(isProfileMode(profile)) {
          if
            constexpr(isProfileWithMiniIntersectionMode(profile)) {
              distinct_intersection_count_ += parent_set_.insert(parent_match).second;
            }
          total_intersection_input_size_ +=
              candidate_set_.size() + ((const G*)current_data_graph_)->getVertexOutDegreeWithHint(parent_match, 0, 0);
          total_intersection_output_size_ += current_targets_.size();
        }
    }
    // intersect(*candidates_, current_data_graph_->getOutNeighbors(parent_match), &current_targets_,
    // input.getKeyMap());
  }
};

class CurrentResults {
 protected:
  const CompressedSubgraphs* input_;
  const uint32_t parent_index_;
  TraverseOperator* owner_;

 public:
  CurrentResults(const CompressedSubgraphs* input, uint32_t parent_index, TraverseOperator* owner)
      : input_(input), parent_index_(parent_index), owner_(owner) {}

  virtual ~CurrentResults() {}

  virtual uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) = 0;
};

template <QueryType profile, typename G>
class CurrentResultsByCandidate : public CurrentResults {
 private:
  uint32_t candidate_index_ = 0;

 public:
  CurrentResultsByCandidate(const CompressedSubgraphs* input, uint32_t parent_index, TraverseOperator* owner)
      : CurrentResults(input, parent_index, owner) {}

  uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    uint32_t n = 0;
    auto& parent_set = input_->getSet(parent_index_);

    auto exceptions = input_->getExceptions(owner_->getSameLabelKeyIndices(), owner_->getSameLabelSetIndices());
    auto data_graph = ((G*)owner_->getCurrentDataGraph());
    for (; n < cap && candidate_index_ < owner_->getCandidateSet()->size(); ++candidate_index_) {
      std::vector<VertexID> parents;
      auto candidate = (*owner_->getCandidateSet())[candidate_index_];
      if (exceptions.count(candidate)) {
        continue;
      }
      intersect(*parent_set, data_graph->getOutNeighborsWithHint(candidate, 0, 0), &parents);  // No need for exceptions
      if
        constexpr(isProfileMode(profile)) {
          owner_->updateIntersectInfo(parent_set->size() + data_graph->getVertexOutDegreeWithHint(candidate, 0, 0),
                                      parents.size());
        }
      if (parents.empty()) {
        continue;
      }
      CompressedSubgraphs output(*input_, parent_index_, makeShared(parents), candidate,
                                 owner_->getSameLabelSetIndices(), owner_->getSetPruningThreshold(), true);
      if (output.empty()) continue;
#ifdef USE_FILTER
      if (owner_->filter(output)) continue;
#endif
      ++n;
      outputs->emplace_back(std::move(output));
    }
    return n;
  }
};

template <QueryType profile, typename G>
class CurrentResultsByParent : public CurrentResults {
  unordered_set<VertexID> exceptions_;

 public:
  CurrentResultsByParent(const CompressedSubgraphs* input, uint32_t parent_index, TraverseOperator* owner)
      : CurrentResults(input, parent_index, owner) {
    exceptions_ = input_->getExceptions(owner_->getSameLabelKeyIndices(), owner_->getSameLabelSetIndices());
  }

  uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    auto& parent_set = *input_->getSet(parent_index_);
    unordered_map<VertexID, int> group_index;
    uint32_t n = 0;
    for (uint32_t i = 0; i < parent_set.size(); ++i) {
      auto parent_match = parent_set[i];
      if (exceptions_.count(parent_match)) continue;
      std::vector<VertexID> targets;
      intersect(*owner_->getCandidateSet(),
                ((G*)owner_->getCurrentDataGraph())->getOutNeighborsWithHint(parent_match, 0, 0), &targets,
                exceptions_);
      if
        constexpr(isProfileMode(profile)) {
          owner_->updateIntersectInfo(
              owner_->getCandidateSet()->size() +
                  ((G*)owner_->getCurrentDataGraph())->getVertexOutDegreeWithHint(parent_match, 0, 0),
              targets.size());
        }
      for (auto target : targets) {
        auto pos = group_index.find(target);
        if (pos == group_index.end()) {
          group_index[target] = outputs->size();
          outputs->emplace_back(*input_, parent_index_, makeVertexSet(parent_match), target,
                                owner_->getSameLabelSetIndices(), owner_->getSetPruningThreshold(), false);
          if (outputs->back().empty()) {
            outputs->pop_back();
            group_index[target] = -1;
          }
          ++n;
        } else if (pos->second != -1) {  // the target is invalid if the group index is -1
          (*outputs)[pos->second].UpdateSet(parent_index_, parent_match);
        }
      }
    }
#ifdef USE_FILTER
    n -= owner_->filter(*outputs, outputs->size() - n, outputs->size());
#endif
    return n;
  }
};

template <QueryType profile, typename G>
class CurrentResultsByExtension : public CurrentResults {
 private:
  unordered_set<VertexID> seen_extensions_;
  std::queue<VertexID> extensions_;
  uint32_t parent_match_index_ = 0;
  unordered_set<VertexID> current_exceptions_;

 public:
  CurrentResultsByExtension(const CompressedSubgraphs* input, uint32_t parent_index, TraverseOperator* owner)
      : CurrentResults(input, parent_index, owner) {
    current_exceptions_ = input_->getExceptions(owner->getSameLabelKeyIndices(), owner->getSameLabelSetIndices());
  }

  uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    uint32_t n = 0;
    auto& parent_set = *input_->getSet(parent_index_);
    while (true) {
      // check existing extensions first
      while (!extensions_.empty()) {
        auto candidate = extensions_.front();
        extensions_.pop();
        std::vector<VertexID> parents;  // valid parents for current candidate
        intersect(parent_set, ((G*)owner_->getCurrentDataGraph())->getOutNeighborsWithHint(candidate, 0, 0),
                  &parents);  // no need for exceptions
        if
          constexpr(isProfileMode(profile)) {
            owner_->updateIntersectInfo(
                parent_set.size() + ((G*)owner_->getCurrentDataGraph())->getVertexOutDegreeWithHint(candidate, 0, 0),
                parents.size());
          }
        CompressedSubgraphs output(*input_, parent_index_, std::make_shared<std::vector<VertexID>>(std::move(parents)),
                                   candidate, owner_->getSameLabelSetIndices(), owner_->getSetPruningThreshold(), true);
        if (output.empty()) {
          continue;
        }
#ifdef USE_FILTER
        if (owner_->filter(output)) continue;
#endif
        outputs->emplace_back(std::move(output));
        if (++n == cap) return n;
      }
      // all parent match extended
      if (parent_match_index_ == parent_set.size()) break;
      // get more extensions by extending from the next parent match
      auto parent_match = parent_set[parent_match_index_];
      getExtensions(parent_match);
      ++parent_match_index_;
    }
    return n;
  }

 private:
  inline void getExtensions(VertexID parent_match) {
    if (current_exceptions_.count(parent_match)) {
      return;
    }
    std::vector<VertexID> current_extensions;
    intersect(*owner_->getCandidateSet(),
              ((G*)owner_->getCurrentDataGraph())->getOutNeighborsWithHint(parent_match, 0, 0), &current_extensions,
              current_exceptions_);
    if
      constexpr(isProfileMode(profile)) {
        owner_->updateIntersectInfo(
            owner_->getCandidateSet()->size() +
                ((G*)owner_->getCurrentDataGraph())->getVertexOutDegreeWithHint(parent_match, 0, 0),
            current_extensions.size());
      }
    for (VertexID neighbor : current_extensions) {
      if (seen_extensions_.insert(neighbor).second) {
        extensions_.push(neighbor);
      }
    }
  }
};

template <typename G>
class ExpandEdgeSetToKeyOperator : public ExpandEdgeOperator {
  enum ExecutionMode { ByCandidate, ByParent, ByExtension };

  uint64_t candidates_neighbor_size_ = 0;
  CurrentResults* current_results_ = nullptr;
  unordered_set<VertexID> parent_set_;  // for profile

 public:
  CONSTRUCT(SetToKey)

  ~ExpandEdgeSetToKeyOperator() { clear(); }

  void clear() {
    if (current_results_ != nullptr) {
      delete current_results_;
      current_results_ = nullptr;
    }
  }

  void input(const std::vector<CompressedSubgraphs>& inputs, const void* data_graph) override {
    // if a different data graph is given, recompute the candidates' neighbor size
    if (data_graph != current_data_graph_) {
      candidates_neighbor_size_ = 0;
      for (auto candidate : *candidates_) {
        candidates_neighbor_size_ += ((G*)data_graph)->getVertexOutDegreeWithHint(candidate, 0, 0);
      }
    }
    TraverseOperator::input(inputs, data_graph);
  }

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    return expandInner<QueryType::Execute>(outputs, cap);
  }

  uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t cap,
                                 uint32_t query_type) override {
    if (query_type == 1) {
      return expandInner<QueryType::Profile>(outputs, cap);
    }
    CHECK_EQ(query_type, 2) << "unknown query type " << query_type;
    return expandInner<QueryType::ProfileWithMiniIntersection>(outputs, cap);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandEdgeSetToKeyOperator";
    toStringInner(ss);
    return ss.str();
  }

  Operator* clone() const override {
    // TODO(tatiana): for now next_ is not handled because it is only used for printing plan
    return new ExpandEdgeSetToKeyOperator(*this);
  }

 private:
  /* Assume the set intersection cost of two sorted sets of size n and m is 2(n+m).
   *
   * cost of enumerating candidates = total set intersection cost
   *                                = 2 * candidates_neighbor_size + 2 * |parent_set| * |candidates_|
   * cost of enumerating parent set
   *     = total set intersection cost + total CompressedSubgraphs lookup cost
   *     = 2 * set_neighbor_size_ + 2 * |parent_set| * |candidates_|
   *       + min(|parent_set| * |candidates_|, set_neighbor_size)
   */
  inline ExecutionMode getExecutionMode(const std::vector<VertexID>* parent_set, uint32_t cap) {
    DCHECK_NE(candidates_neighbor_size_, 0);
    uint64_t set_neighbor_size = 0;
    for (auto v : *parent_set) {
      set_neighbor_size += ((G*)current_data_graph_)->getVertexOutDegreeWithHint(v, 0, 0);
    }
    auto enumerating_candidate_cost = 2 * candidates_neighbor_size_;
    auto enumerating_parent_cost =
        2 * set_neighbor_size + std::min(parent_set->size() * candidates_->size(), set_neighbor_size);
    if (enumerating_candidate_cost < enumerating_parent_cost) return ByCandidate;

    /* ByParent requires all parent to be processed before the outputs become ready, and thus the output size is unknown
     * in advance and unbounded. ByExtension incurs extra set intersection, but the output size can be bounded. */
    if (std::min(candidates_->size(), set_neighbor_size) <=
        cap) {  // here we use set_neighbor_size to approximate the output size. as long as the
                // output size does not exceed cap, ByParent is preferred than ByExtension.
                // TODO(tatiana): consider better estimation on the output size
      return ByParent;
    }
    return ByExtension;
  }

  template <QueryType profile>
  uint32_t expandInner(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) {
    DCHECK_GT(cap, 0);
    uint32_t needed = cap;
    while (true) {
      if (current_results_ != nullptr) {
        uint32_t got = current_results_->getResults(outputs, needed);
        DCHECK_LE(got, needed);
        needed -= got;

        if (needed == 0) {
          return cap - needed;
        }

        if
          constexpr(isProfileWithMiniIntersectionMode(profile)) {
            auto& parent_set = *(*current_inputs_)[input_index_ - 1].getSet(parent_index_);
            for (auto parent_match : parent_set) {
              distinct_intersection_count_ += parent_set_.insert(parent_match).second;
            }
          }
        delete current_results_;
        current_results_ = nullptr;
      }
      if (input_index_ == current_inputs_->size()) {
        return cap - needed;
      }
      expandInner<profile>((*current_inputs_)[input_index_], needed);
      if
        constexpr(isProfileMode(profile)) {
          total_num_input_subgraphs_ += (*current_inputs_)[input_index_].getNumSubgraphs();
        }
      ++input_index_;
    }
    return cap - needed;
  }

  template <QueryType profile>
  void expandInner(const CompressedSubgraphs& input, uint32_t cap) {
    auto& parent_set = input.getSet(parent_index_);
    DCHECK(current_results_ == nullptr);
    ExecutionMode mode = getExecutionMode(parent_set.get(), cap);
    switch (mode) {
    case ByCandidate:
      // enumerating target candidates is likely to incur less computation
      current_results_ = new CurrentResultsByCandidate<profile, G>(&input, parent_index_, this);
      break;
    case ByParent:
      // enumerating vertices in the parent set is likely to incur less computation
      current_results_ = new CurrentResultsByParent<profile, G>(&input, parent_index_, this);
      break;
    default:
      // enumerating vertices in the parent set is likely to incur less computation
      current_results_ = new CurrentResultsByExtension<profile, G>(&input, parent_index_, this);
    }
  }
};

TraverseOperator* ExpandEdgeOperator::newExpandEdgeKeyToSetOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex, const unordered_map<QueryVertexID, uint32_t>& indices,
    const std::vector<uint32_t>& same_label_key_indices, const std::vector<uint32_t>& same_label_set_indices,
    uint64_t set_pruning_threshold, SubgraphFilter* filter, GraphType graph_type) {
  DCHECK_GT(indices.count(parent_vertex), 0);
  DCHECK_GT(indices.count(target_vertex), 0);
  if (graph_type == GraphType::Normal) {
    return new ExpandEdgeKeyToSetOperator<Graph>(indices.at(parent_vertex), indices.at(target_vertex), parent_vertex,
                                                 target_vertex, same_label_key_indices, same_label_set_indices,
                                                 set_pruning_threshold, filter);
  }
  if (graph_type == GraphType::GraphView) {
    return new ExpandEdgeKeyToSetOperator<GraphView<Graph>>(indices.at(parent_vertex), indices.at(target_vertex),
                                                            parent_vertex, target_vertex, same_label_key_indices,
                                                            same_label_set_indices, set_pruning_threshold, filter);
  }
  CHECK(graph_type == GraphType::BipartiteGraphView) << "unknown graph type " << ((uint32_t)graph_type);
  return new ExpandEdgeKeyToSetOperator<GraphView<BipartiteGraph>>(
      indices.at(parent_vertex), indices.at(target_vertex), parent_vertex, target_vertex, same_label_key_indices,
      same_label_set_indices, set_pruning_threshold, filter);
}

TraverseOperator* ExpandEdgeOperator::newExpandEdgeKeyToKeyOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex, const unordered_map<QueryVertexID, uint32_t>& indices,
    const std::vector<uint32_t>& same_label_key_indices, const std::vector<uint32_t>& same_label_set_indices,
    uint64_t set_pruning_threshold, SubgraphFilter* filter, bool intersect_candidates, GraphType graph_type) {
  DCHECK_GT(indices.count(parent_vertex), 0);
  DCHECK_GT(indices.count(target_vertex), 0);
  if (graph_type == GraphType::Normal) {
    if (intersect_candidates) {
      return new ExpandEdgeKeyToKeyOperator<Graph, true>(indices.at(parent_vertex), indices.at(target_vertex),
                                                         parent_vertex, target_vertex, same_label_key_indices,
                                                         same_label_set_indices, set_pruning_threshold, filter);
    }
    return new ExpandEdgeKeyToKeyOperator<Graph, false>(indices.at(parent_vertex), indices.at(target_vertex),
                                                        parent_vertex, target_vertex, same_label_key_indices,
                                                        same_label_set_indices, set_pruning_threshold, filter);
  }
  if (graph_type == GraphType::GraphView) {
    if (intersect_candidates) {
      return new ExpandEdgeKeyToKeyOperator<GraphView<Graph>, true>(
          indices.at(parent_vertex), indices.at(target_vertex), parent_vertex, target_vertex, same_label_key_indices,
          same_label_set_indices, set_pruning_threshold, filter);
    }
    return new ExpandEdgeKeyToKeyOperator<GraphView<Graph>, false>(
        indices.at(parent_vertex), indices.at(target_vertex), parent_vertex, target_vertex, same_label_key_indices,
        same_label_set_indices, set_pruning_threshold, filter);
  }
  CHECK(graph_type == GraphType::BipartiteGraphView) << "unknown graph type " << ((uint32_t)graph_type);
  {
    if (intersect_candidates) {
      return new ExpandEdgeKeyToKeyOperator<GraphView<BipartiteGraph>, true>(
          indices.at(parent_vertex), indices.at(target_vertex), parent_vertex, target_vertex, same_label_key_indices,
          same_label_set_indices, set_pruning_threshold, filter);
    }
    return new ExpandEdgeKeyToKeyOperator<GraphView<BipartiteGraph>, false>(
        indices.at(parent_vertex), indices.at(target_vertex), parent_vertex, target_vertex, same_label_key_indices,
        same_label_set_indices, set_pruning_threshold, filter);
  }
}

// tricky case: expand, look up group for each match of target, and copy parent to set for each group; or
// consider an alternative: enumerate each candidate of target, and do set intersection between the target
// neighbors and parent sets
TraverseOperator* ExpandEdgeOperator::newExpandEdgeSetToKeyOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex, const unordered_map<QueryVertexID, uint32_t>& indices,
    const std::vector<uint32_t>& same_label_key_indices, const std::vector<uint32_t>& same_label_set_indices,
    uint64_t set_pruning_threshold, SubgraphFilter* filter, GraphType graph_type) {
  DCHECK_GT(indices.count(parent_vertex), 0);
  DCHECK_GT(indices.count(target_vertex), 0);
  if (graph_type == GraphType::Normal) {
    return new ExpandEdgeSetToKeyOperator<Graph>(indices.at(parent_vertex), indices.at(target_vertex), parent_vertex,
                                                 target_vertex, same_label_key_indices, same_label_set_indices,
                                                 set_pruning_threshold, filter);
  }
  if (graph_type == GraphType::GraphView) {
    return new ExpandEdgeSetToKeyOperator<GraphView<Graph>>(indices.at(parent_vertex), indices.at(target_vertex),
                                                            parent_vertex, target_vertex, same_label_key_indices,
                                                            same_label_set_indices, set_pruning_threshold, filter);
  }
  CHECK(graph_type == GraphType::BipartiteGraphView) << "unknown graph type " << ((uint32_t)graph_type);
  return new ExpandEdgeSetToKeyOperator<GraphView<BipartiteGraph>>(
      indices.at(parent_vertex), indices.at(target_vertex), parent_vertex, target_vertex, same_label_key_indices,
      same_label_set_indices, set_pruning_threshold, filter);
}

}  // namespace circinus

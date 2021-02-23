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
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/types.h"
#include "ops/traverse_operator.h"

namespace circinus {

#define makeVertexSet(vertex) std::make_shared<std::vector<VertexID>>(std::vector<VertexID>({vertex})

class ExpandEdgeKeyToSetOperator : public TraverseOperator {
  uint32_t parent_index_;  // index of parent query vertex in the compressed subgraphs
  uint32_t target_index_;  // index of target query vertex in the compressed subgraphs

  // TODO(tatiana): statistics for how much saving is made

 public:
  ExpandEdgeKeyToSetOperator(uint32_t parent_index, uint32_t target_index)
      : parent_index_(parent_index), target_index_(target_index) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    uint32_t n = 0;
    for (; n < cap && input_index_ < current_inputs_->size(); ++input_index_) {
      n += expandInner(outputs, (*current_inputs_)[input_index_]);
    }
    return n;
  }

 private:
  /** @returns True if one CompressedSubgraphs is generated, else false. */
  inline bool expandInner(std::vector<CompressedSubgraphs>* outputs, const CompressedSubgraphs& input) {
    std::vector<VertexID> targets;
    auto parent_match = input.getKeyVal(parent_index_);
    intersect(*candidates_, current_data_graph_->getOutNeighbors(parent_match), &targets);
    if (targets.empty()) {
      return false;
    }
    outputs->emplace_back(input, std::make_shared<std::vector<VertexID>>(std::move(targets)));
    return true;
  }
};

class ExpandEdgeKeyToKeyOperator : public TraverseOperator {
  uint32_t parent_index_;  // index of parent query vertex in the compressed subgraphs
  uint32_t target_index_;  // index of target query vertex in the compressed subgraphs

  // calculated from current_inputs_[input_index_]
  std::vector<VertexID> current_targets_;
  uint32_t current_target_index_ = 0;

 public:
  ExpandEdgeKeyToKeyOperator(uint32_t parent_index, uint32_t target_index)
      : parent_index_(parent_index), target_index_(target_index) {}

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    uint32_t n = 0;
    while (true) {
      // if there are existing targets from the last input, consume first
      if (current_target_index_ < current_targets_.size()) {
        auto target_end = std::min(current_target_index_ + cap - n, (uint32_t)current_targets_.size());
        n += target_end - current_target_index_;
        for (; current_target_index_ < target_end; ++current_target_index_) {
          outputs->emplace_back((*current_inputs_)[input_index_ - 1], current_targets_[current_target_index_]);
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
      expandInner((*current_inputs_)[input_index_]);
      ++input_index_;
    }
    return n;
  }

 private:
  inline void expandInner(const CompressedSubgraphs& input) {
    current_targets_.clear();
    current_target_index_ = 0;
    auto parent_match = input.getKeyVal(parent_index_);
    intersect(*candidates_, current_data_graph_->getOutNeighbors(parent_match), &current_targets_);
  }
};

class CurrentResults {
 protected:
  const std::vector<VertexID>* candidates_;
  const CompressedSubgraphs* input_;
  const Graph* data_graph_;
  const uint32_t parent_index_;

 public:
  CurrentResults(const std::vector<VertexID>* candidates, const CompressedSubgraphs* input, const Graph* data_graph,
                 uint32_t parent_index)
      : candidates_(candidates), input_(input), data_graph_(data_graph), parent_index_(parent_index) {}

  virtual ~CurrentResults() {}

  virtual uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) = 0;
};

class CurrentResultsByCandidate : public CurrentResults {
 private:
  uint32_t candidate_index_ = 0;

 public:
  CurrentResultsByCandidate(const std::vector<VertexID>* candidates, const CompressedSubgraphs* input,
                            const Graph* data_graph, uint32_t parent_index)
      : CurrentResults(candidates, input, data_graph, parent_index) {}

  uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    uint32_t n = 0;
    auto& parent_set = input_->getSet(parent_index_);
    for (; n < cap && candidate_index_ < candidates_->size(); ++candidate_index_) {
      std::vector<VertexID> parents;
      auto candidate = (*candidates_)[candidate_index_];
      intersect(*parent_set, data_graph_->getOutNeighbors(candidate), &parents);
      if (parents.empty()) {
        continue;
      }
      ++n;
      outputs->emplace_back(*input_, parent_index_, std::make_shared<std::vector<VertexID>>(std::move(parents)),
                            candidate);
    }
    return n;
  }
};

class CurrentResultsByParent : public CurrentResults {
 public:
  CurrentResultsByParent(const std::vector<VertexID>* candidates, const CompressedSubgraphs* input,
                         const Graph* data_graph, uint32_t parent_index)
      : CurrentResults(candidates, input, data_graph, parent_index) {}

  uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    auto& parent_set = *input_->getSet(parent_index_);
    std::unordered_map<VertexID, uint32_t> group_index;
    uint32_t n = 0;
    for (uint32_t i = 0; i < parent_set.size(); ++i) {
      auto parent_match = parent_set[i];
      std::vector<VertexID> targets;
      intersect(*candidates_, data_graph_->getOutNeighbors(parent_match), &targets);
      for (auto target : targets) {
        auto pos = group_index.find(target);
        if (pos == group_index.end()) {
          group_index[target] = outputs->size();
          outputs->emplace_back(*input_, parent_index_, makeVertexSet(parent_match)), target);
          ++n;
        } else {
          (*outputs)[pos->second].UpdateSet(parent_index_, parent_match);
        }
      }
    }
    return n;
  }
};

class CurrentResultsByExtension : public CurrentResults {
 private:
  std::unordered_set<VertexID> seen_extensions_;
  std::vector<VertexID> extensions_;
  uint32_t parent_match_index_ = 0;

 public:
  CurrentResultsByExtension(const std::vector<VertexID>* candidates, const CompressedSubgraphs* input,
                            const Graph* data_graph, uint32_t parent_index)
      : CurrentResults(candidates, input, data_graph, parent_index) {}

  uint32_t getResults(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    getExtensions(cap);

    auto& parent_set = *input_->getSet(parent_index_);
    uint32_t n = std::min(cap, (uint32_t)extensions_.size());
    int leftover_size = extensions_.size() - n;
    for (int i = extensions_.size() - 1; i >= leftover_size; --i) {
      auto candidate = extensions_[i];
      std::vector<VertexID> parents;  // valid parents for current candidate
      intersect(parent_set, data_graph_->getOutNeighbors(candidate), &parents);
      outputs->emplace_back(*input_, parent_index_, std::make_shared<std::vector<VertexID>>(std::move(parents)),
                            candidate);
    }
    extensions_.resize(leftover_size);
    return n;
  }

 private:
  inline void getExtensions(uint32_t cap) {
    auto& parent_set = *input_->getSet(parent_index_);
    std::vector<VertexID> current_extensions;
    for (; extensions_.size() < cap && parent_match_index_ < parent_set.size(); ++parent_match_index_) {
      auto parent_match = parent_set[parent_match_index_];
      auto neighbors = data_graph_->getOutNeighbors(parent_match);
      current_extensions.reserve(neighbors.second);
      for (uint32_t i = 0; i < neighbors.second; ++i) {
        auto neighbor = neighbors.first[i];
        if (seen_extensions_.insert(neighbor).second) {
          current_extensions.push_back(neighbor);
        }
      }
      intersect(*candidates_, current_extensions, &extensions_);
      current_extensions.clear();
    }
  }
};

class ExpandEdgeSetToKeyOperator : public TraverseOperator {
  enum ExecutionMode { ByCandidate, ByParent, ByExtension };

  uint32_t parent_index_;  // index of parent query vertex in the compressed subgraphs
  uint32_t target_index_;  // index of target query vertex in the compressed subgraphs

  uint64_t candidates_neighbor_size_ = 0;

  CurrentResults* current_results_ = nullptr;

 public:
  ExpandEdgeSetToKeyOperator(uint32_t parent_index, uint32_t target_index)
      : parent_index_(parent_index), target_index_(target_index) {}

  ~ExpandEdgeSetToKeyOperator() { clear(); }

  void clear() {
    if (current_results_ != nullptr) {
      delete current_results_;
      current_results_ = nullptr;
    }
  }

  void input(const std::vector<CompressedSubgraphs>& inputs, const Graph* data_graph) override {
    // if a different data graph is given, recompute the candidates' neighbor size
    if (data_graph != current_data_graph_) {
      candidates_neighbor_size_ = 0;
      for (auto candidate : *candidates_) {
        candidates_neighbor_size_ += data_graph->getVertexOutDegree(candidate);
      }
    }
    TraverseOperator::input(inputs, data_graph);
  }

  uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) override {
    DCHECK_GT(cap, 0);
    uint32_t needed = cap;
    while (true) {
      if (current_results_ != nullptr) {
        auto got = current_results_->getResults(outputs, needed);
        DCHECK_LE(got, needed);
        needed -= got;
        if (needed == 0) {
          return cap;
        }
        delete current_results_;
        current_results_ = nullptr;
      }
      if (input_index_ == current_inputs_->size()) {
        return cap - needed;
      }
      expandInner((*current_inputs_)[input_index_], needed);
      ++input_index_;
    }
    return cap - needed;
  }

 private:
  /* Assume the set intersection cost of two sorted sets of size n and m is 2(n+m).
   *
   * cost of enumerating candidates = total set intersection cost
   *                                = 2 * set_neighbor_size + 2 * |parent_set| * |candidates_|
   * cost of enumerating parent set
   *     = total set intersection cost + total CompressedSubgraphs lookup cost
   *     = 2 * candidates_neighbor_size_ + 2 * |parent_set| * |candidates_|
   *       + min(|parent_set| * |candidates_|, set_neighbor_size)
   */
  inline ExecutionMode getExecutionMode(const std::vector<VertexID>* parent_set, uint32_t cap) {
    uint64_t set_neighbor_size = 0;
    for (auto v : *parent_set) {
      set_neighbor_size += current_data_graph_->getVertexOutDegree(v);
    }
    auto enumerating_candidate_cost = 2 * candidates_neighbor_size_;
    auto enumerating_parent_cost =
        2 * set_neighbor_size + std::min(parent_set->size() * candidates_->size(), set_neighbor_size);
    if (enumerating_candidate_cost < enumerating_parent_cost) return ByCandidate;

    /* ByParent requires all parent to be processed before the outputs become ready, and thus the output size is unknown
     * in advance and unbounded. ByExtension incurs extra set intersection, but the output size can be bounded. */
    if (set_neighbor_size <= cap) {  // here we use set_neighbor_size to approximate the output size. as long as the
                                     // output size does not exceed cap, ByParent is preferred than ByExtension.
                                     // TODO(tatiana): consider better estimation on the output size
      return ByParent;
    }
    return ByExtension;
  }

  void expandInner(const CompressedSubgraphs& input, uint32_t cap) {
    auto& parent_set = input.getSet(parent_index_);
    DCHECK(current_results_ == nullptr);
    auto mode = getExecutionMode(parent_set.get(), cap);
    switch (mode) {
    case ByCandidate:
      // enumerating target candidates is likely to incur less computation
      current_results_ = new CurrentResultsByCandidate(candidates_, &input, current_data_graph_, parent_index_);
      break;
    case ByParent:
      // enumerating vertices in the parent set is likely to incur less computation
      current_results_ = new CurrentResultsByParent(candidates_, &input, current_data_graph_, parent_index_);
      break;
    default:
      // enumerating vertices in the parent set is likely to incur less computation
      current_results_ = new CurrentResultsByExtension(candidates_, &input, current_data_graph_, parent_index_);
    }
  }
};

TraverseOperator* ExpandEdgeOperator::newExpandEdgeOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex, const std::vector<int>& cover_table,
    const std::unordered_map<QueryVertexID, uint32_t>& indices) {
  CHECK_GT(indices.count(parent_vertex), 0);
  CHECK_GT(indices.count(target_vertex), 0);
  // the target is not a compression key, and the parent must be in the cover: expand and copy target list
  if (cover_table[target_vertex] != 1) {
    return new ExpandEdgeKeyToSetOperator(indices.at(parent_vertex), indices.at(target_vertex));
  }
  // the target is a compression key
  if (cover_table[parent_vertex] == 1) {
    // the parent is a compression key: expand and enumerate parent-target pairs
    return new ExpandEdgeKeyToKeyOperator(indices.at(parent_vertex), indices.at(target_vertex));
  }

  // tricky case: expand, look up group for each match of target, and copy parent to set for each group; or
  // consider an alternative: enumerate each candidate of target, and do set intersection between the target
  // neighbors and parent sets
  return new ExpandEdgeSetToKeyOperator(indices.at(parent_vertex), indices.at(target_vertex));
}

}  // namespace circinus

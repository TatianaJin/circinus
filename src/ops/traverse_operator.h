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
#include <chrono>
#include <string>
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/bipartite_graph.h"
#include "graph/types.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/operator.h"
#include "utils/hashmap.h"
#include "utils/utils.h"

namespace circinus {

inline void removeExceptions(const std::pair<const VertexID*, uint32_t>& setPair, std::vector<VertexID>* res,
                             const unordered_set<VertexID>& except = {}) {
  for(uint32_t i = 0; i < setPair.second; ++i)
  {
    auto vid = setPair.first[i];
    if(except.count(vid) == 0) res->emplace_back(vid);
  }
}


/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::pair<const VertexID*, uint32_t>& set1,
                      const std::pair<const VertexID*, uint32_t>& set2, std::vector<VertexID>* intersection,
                      const unordered_set<VertexID>& except = {}) {
  if (set1.second <= set2.second) {
    intersection->reserve(set1.second);
    uint32_t set2_index = 0;
    for (uint32_t i = 0; i < set1.second; ++i) {
      auto vid = set1.first[i];
      while (set2_index < set2.second && set2.first[set2_index] < vid) {
        ++set2_index;
      }
      if (set2_index == set2.second) {  // all elements in the rest of set2 are smaller than the rest of set1
        break;
      }
      if (set2.first[set2_index] == vid && except.count(vid) == 0) intersection->emplace_back(vid);
    }
  } else {
    intersect(set2, set1, intersection, except);
  }
}

inline void intersect(const unordered_set<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection, const unordered_set<VertexID>& except = {}) {
  intersection->reserve(std::min(set2.second, (uint32_t)set1.size()));
  for (uint32_t i = 0; i < set2.second; ++i) {
    auto vid = set2.first[i];
    if (set1.count(vid) && except.count(vid) == 0) intersection->emplace_back(vid);
  }
}

/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection, const unordered_set<VertexID>& except = {}) {
  intersect(std::make_pair(set1.data(), (uint32_t)set1.size()), set2, intersection, except);
}

/** set1 and set2 must be sorted in ascending order */
inline void intersect(const std::vector<VertexID>& set1, const std::vector<VertexID>& set2,
                      std::vector<VertexID>* intersection, const unordered_set<VertexID>& except = {}) {
  intersect(set1, std::make_pair(set2.data(), (uint32_t)set2.size()), intersection, except);
}

void intersectInplace(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection);

class TraverseOperator : public Operator {
 protected:
  const std::vector<VertexID>* candidates_ = nullptr;

  /* for non-repeated-vertex check */
  uint64_t set_pruning_threshold_ = ~0u;
  std::vector<uint32_t> same_label_key_indices_;
  std::vector<uint32_t> same_label_set_indices_;
  SubgraphFilter* const subgraph_filter_ = nullptr;  // owned by the execution plan

  /* transient variables for recording the current inputs */
  uint32_t input_index_ = 0;
  const std::vector<CompressedSubgraphs>* current_inputs_ = nullptr;
  const Graph* current_data_graph_ = nullptr;

  uint32_t current_inputs_size_ = 0;
  /* for profiling */
  uint64_t total_input_size_ = 0;
  uint64_t total_output_size_ = 0;
  uint64_t total_num_input_subgraphs_ = 0;
  uint64_t total_num_output_subgraphs_ = 0;
  double total_time_in_milliseconds_ = 0;
  // to be updated in derived class
  uint64_t intersection_count_ = 0;
  uint64_t total_intersection_input_size_ = 0;
  uint64_t total_intersection_output_size_ = 0;
  uint64_t distinct_intersection_count_ =
      0;  // the minimal number of intersection needed if all intersection function call results can be cached
  std::vector<BipartiteGraph*> bg_pointers_;
  bool use_bipartite_graph_flag=false;

 public:
  TraverseOperator() {}
  explicit TraverseOperator(SubgraphFilter* filter) : subgraph_filter_(filter) {}
  TraverseOperator(const std::vector<uint32_t>& same_label_key_indices,
                   const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,
                   SubgraphFilter* filter = nullptr)
      : set_pruning_threshold_(set_pruning_threshold),
        same_label_key_indices_(same_label_key_indices),
        same_label_set_indices_(same_label_set_indices),
        subgraph_filter_(filter) {}
  virtual ~TraverseOperator() {}

  inline virtual void setCandidateSets(const std::vector<VertexID>* candidates) { candidates_ = candidates; }
  inline const std::vector<VertexID>* getCandidateSets() const { return candidates_; }
  inline const uint32_t getInputIndex() const { return input_index_; }
  inline const auto& getSameLabelKeyIndices() const { return same_label_key_indices_; }
  inline const auto& getSameLabelSetIndices() const { return same_label_set_indices_; }
  inline auto getSetPruningThreshold() const { return set_pruning_threshold_; }

  inline bool filter(const CompressedSubgraphs& subgraphs) {
    DCHECK(subgraph_filter_ != nullptr);
    return subgraph_filter_->filter(subgraphs);
  }

  inline bool filter(std::vector<CompressedSubgraphs>& subgraphs, uint32_t start, uint32_t end) {
    DCHECK(subgraph_filter_ != nullptr);
    return subgraph_filter_->filter(subgraphs, start, end);
  }

  virtual void input(const std::vector<CompressedSubgraphs>& inputs, const Graph* data_graph) {
    current_inputs_ = &inputs;
    input_index_ = 0;
    current_data_graph_ = data_graph;
  }

  virtual uint32_t expand(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) = 0;

  virtual void inputAndProfile(const std::vector<CompressedSubgraphs>& inputs, const Graph* data_graph) {
    input(inputs, data_graph);
    current_inputs_size_ = inputs.size();
    total_input_size_ += inputs.size();
  }

  uint32_t expandAndProfile(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) {
    auto start = std::chrono::high_resolution_clock::now();
    auto n = expandAndProfileInner(outputs, cap);
    auto stop = std::chrono::high_resolution_clock::now();
    total_time_in_milliseconds_ +=
        (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000000.0);
    {
      uint32_t size = 0;
      auto offset = outputs->size() - n;
      for (uint32_t i = 0; i < n; ++i) {
        if ((*outputs)[offset + i].getNumIsomorphicSubgraphs(1) == 0) {
          CHECK(false) << toString() << "\n\t" << (*outputs)[offset + i].toString();
          continue;
        }
        if (size != i) {
          (*outputs)[offset + size] = std::move((*outputs)[offset + i]);
        }
        ++size;
      }
      outputs->erase(outputs->begin() + (offset + size), outputs->end());
      n = size;
    }
    total_num_output_subgraphs_ += getNumSubgraphs(*outputs, outputs->size() - n, outputs->size());
    total_output_size_ += n;
    return n;
  }

  void updateIntersectInfo(uint32_t input_size, uint32_t output_size) {
    ++intersection_count_;
    total_intersection_input_size_ += input_size;
    total_intersection_output_size_ += output_size;
  }

  std::string toProfileString() const override {
    std::stringstream ss;
    ss << toString() << ',' << total_time_in_milliseconds_ << ',' << getTotalInputSize() << ',' << total_output_size_
       << ',' << total_num_input_subgraphs_ << ',' << total_num_output_subgraphs_ << ',' << intersection_count_ << ','
       << total_intersection_input_size_ << ',' << total_intersection_output_size_ << ','
       << distinct_intersection_count_;
    return ss.str();
  }

  void addBipartiteGraph(BipartiteGraph* bg)
  {
    bg_pointers_.emplace_back(bg);
  }

  void useBipartiteGraph(std::vector<std::vector<VertexID>> candidate_sets) // must used after input() called
  {
    use_bipartite_graph_flag=1;
    for(auto p:bg_pointers_)
    {
      p->populateGraph(current_data_graph_,candidate_sets);
    }
  }

  void intersectOrNot(const std::vector<VertexID>& set1, const std::pair<const VertexID*, uint32_t>& set2,
                      std::vector<VertexID>* intersection, const unordered_set<VertexID>& except = {}) {
    if(use_bipartite_graph_flag)removeExceptions(set2, intersection,except);
    else intersect(set1, set2, intersection,except);
  }

 protected:
  virtual uint32_t expandAndProfileInner(std::vector<CompressedSubgraphs>* outputs, uint32_t cap) = 0;
  inline uint32_t getTotalInputSize() const { return total_input_size_ - (current_inputs_size_ - input_index_); }
};

}  // namespace circinus

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

#include "ops/output_operator.h"

#include <fstream>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "glog/logging.h"

#include "graph/compressed_subgraphs.h"
#include "ops/operator.h"

namespace circinus {

class CountOutputOperator : public OutputOperator {
 public:
  explicit CountOutputOperator(SameLabelIndices&& same_label_indices) : OutputOperator(std::move(same_label_indices)) {}

  bool validateAndOutput(const std::vector<CompressedSubgraphs>& input, uint32_t& input_start, uint32_t input_end,
                         uint32_t output_index) const override {
    auto count_acc = outputs_->getCount(output_index);
    if (!set_less_than_constraints_.empty()) {
      while (input_start < input_end) {
        auto& group = input[input_start];
        auto update = group.getNumIsomorphicSubgraphsWithConstraints(
            same_label_indices_, constraints_adjs_, enumerate_orders_, outputs_->getLimitPerThread() - count_acc);
        count_acc = outputs_->updateCount(update, output_index);
        ++input_start;
        if (count_acc >= outputs_->getLimitPerThread()) {
          return true;
        }
      }
      return false;
    }
    while (input_start < input_end) {
      auto& group = input[input_start];
      auto update = group.getNumIsomorphicSubgraphs(same_label_indices_, outputs_->getLimitPerThread() - count_acc);
      count_acc = outputs_->updateCount(update, output_index);
      ++input_start;
      if (count_acc >= outputs_->getLimitPerThread()) {
        // DLOG(INFO) << "last input num subgraphs " << group.getNumSubgraphs() << " isomorphic " << update << " total "
        //            << count_acc << " limit " << outputs_->getLimitPerThread();
        return true;
      }
    }
    return false;
  }

  std::string toString() const override { return "CountOutputOperator"; }
};

void OutputOperator::setPartialOrder(std::vector<std::pair<uint32_t, uint32_t>>&& constraints) {
  if (constraints.empty()) return;
  LOG(INFO) << "start set partial order";
  auto set_group_size = same_label_indices_.size();
  enumerate_orders_.resize(set_group_size);
  constraints_adjs_.resize(set_group_size);
  for (uint32_t i = 0; i < set_group_size; ++i) {
    auto& same_label_set_index = same_label_indices_[i].second;
    uint32_t set_index_size = same_label_set_index.size();

    // treat each a<b constraint as an edge a->b
    unordered_map<uint32_t, uint32_t> same_label_index_set;  // set index to same-label group index
    std::vector<uint32_t> in_degree(set_index_size, 0);
    std::vector<std::vector<uint32_t>> tmp_adj(set_index_size);
    std::vector<uint32_t> index_to_order(set_index_size);  // mapping group index to enumerate order index
    constraints_adjs_[i].resize(set_index_size);
    for (uint32_t j = 0; j < set_index_size; ++j) {
      same_label_index_set.insert({same_label_set_index[j], j});
    }

    // construct constraint graph
    for (auto constraint : constraints) {
      auto pos = same_label_index_set.find(constraint.first);
      if (pos != same_label_index_set.end()) {
        auto dst = same_label_index_set[constraint.second];
        ++in_degree[dst];
        tmp_adj[pos->second].push_back(dst);
      }
    }

    // get enumerate order by topological order following the edge direction
    std::queue<uint32_t> current_top_vertices;
    for (uint32_t j = 0; j < set_index_size; ++j) {
      if (in_degree[j] == 0) {
        current_top_vertices.push(j);
      }
    }
    while (!current_top_vertices.empty()) {
      uint32_t cur_vertex = current_top_vertices.front();
      current_top_vertices.pop();
      index_to_order[cur_vertex] = enumerate_orders_[i].size();
      enumerate_orders_[i].push_back(same_label_set_index[cur_vertex]);
      for (uint32_t out_nbr : tmp_adj[cur_vertex]) {
        if (--in_degree[out_nbr] == 0) {
          current_top_vertices.push(out_nbr);
        }
      }
    }
    CHECK_EQ(enumerate_orders_[i].size(), set_index_size);

    // get constraint adj, have mapped the index to order
    for (uint32_t j = 0; j < set_index_size; ++j) {
      auto tmp_src = same_label_index_set.at(enumerate_orders_[i][j]);
      for (auto tmp_dst : tmp_adj[tmp_src]) {
        constraints_adjs_[i][index_to_order[tmp_dst]].push_back(j);
      }
    }
  }
  set_less_than_constraints_ = std::move(constraints);
  LOG(INFO) << "set less than constraints, " << set_less_than_constraints_.size();
  LOG(INFO) << "end set partial order";
}

Outputs& Outputs::init(uint32_t n_threads, const std::string& path_prefix) {
  n_threads_ = n_threads;
  if (path_prefix != "") {
    output_file_per_thread_.resize(n_threads);
    for (uint32_t i = 0; i < n_threads; ++i) {
      output_file_per_thread_[i].open(path_prefix + '_' + std::to_string(i));
      CHECK(output_file_per_thread_[i].is_open()) << "Cannot write to " << path_prefix << '_' << i;
    }
  } else {
    n_matches_per_thread_.resize(n_threads, 0);
  }
  return *this;
}

OutputOperator* OutputOperator::newOutputOperator(OutputType type, SameLabelIndices&& same_label_indices) {
  CHECK(type == OutputType::Count) << "unsupported output type " << ((uint32_t)type);
  return new CountOutputOperator(std::move(same_label_indices));
}

}  // namespace circinus

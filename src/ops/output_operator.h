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

#include <chrono>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/compressed_subgraphs.h"
#include "ops/operator.h"
#include "utils/utils.h"

namespace circinus {

enum class OutputType : uint32_t { Subgraph = 0, Count }; /* TODO(tatiana): support enumeration */

class Outputs {
  std::vector<uint64_t> n_matches_per_thread_;
  std::vector<std::ofstream> output_file_per_thread_;
  uint32_t n_threads_ = 0;
  uint64_t limit_per_thread_ = ~0u;

 public:
  Outputs& init(uint32_t n_threads, const std::string& path_prefix = "");

  inline uint64_t updateCount(uint64_t count, uint32_t thread_id) {
    DCHECK_LT(thread_id, n_matches_per_thread_.size());
    n_matches_per_thread_[thread_id] += count;
    return n_matches_per_thread_[thread_id];
  }

  inline uint64_t getCount(uint32_t thread_id) const { return n_matches_per_thread_[thread_id]; }

  // FIXME(tatiana): problematic when running multiple threads, and the outputs are all in the task running on one
  // thread but get skipped due to the limit
  inline void limit(uint64_t total_limit) {
    limit_per_thread_ = total_limit / n_threads_ + ((total_limit % n_threads_) != 0);
  }

  inline uint64_t getLimitPerThread() const { return limit_per_thread_; }

  inline uint64_t getCount() const {
    uint64_t total = 0;
    for (auto count : n_matches_per_thread_) {
      total += count;
    }
    return total;
  }
};

class OutputOperator : public Operator {
 protected:
  using SameLabelIndices = std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>;  // {{keys}, {sets}}

  Outputs* outputs_ = nullptr;
  SameLabelIndices same_label_indices_;

  double total_time_in_milliseconds_ = 0;
  uint64_t total_input_size_ = 0;
  uint64_t total_num_input_subgraphs_ = 0;
  uint64_t leftover_input_ = 0;

 public:
  explicit OutputOperator(SameLabelIndices&& same_label_indices) : same_label_indices_(std::move(same_label_indices)) {}
  virtual ~OutputOperator() {}

  /**
   * @param same_label_indices The indices of vertices of the same label {{keys}, {sets}}
   */
  static OutputOperator* newOutputOperator(
      OutputType type, std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>&& same_label_indices);

  void setOutput(Outputs* outputs) { outputs_ = outputs; }

  virtual bool validateAndOutput(const std::vector<CompressedSubgraphs>& input, uint32_t output_index) = 0;

  inline bool validateAndOutputAndProfile(const std::vector<CompressedSubgraphs>& input, uint32_t output_index) {
    auto start = std::chrono::high_resolution_clock::now();
    auto ret = validateAndOutput(input, output_index);
    auto stop = std::chrono::high_resolution_clock::now();
    total_time_in_milliseconds_ +=
        (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000000.0);
    total_num_input_subgraphs_ += getNumSubgraphs(input, 0, input.size() - leftover_input_);
    total_input_size_ += input.size() - leftover_input_;
    return ret;
  }

  std::string toProfileString() const override {
    std::stringstream ss;
    ss << toString() << ',' << total_time_in_milliseconds_ << ',' << total_input_size_ << ",,"
       << total_num_input_subgraphs_;
    return ss.str();
  }
};

}  // namespace circinus

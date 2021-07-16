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
#include <string>
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

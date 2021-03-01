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

#include <string>
#include <vector>

#include "glog/logging.h"

#include "graph/compressed_subgraphs.h"
#include "ops/operator.h"

namespace circinus {

class OutputOperator : public Operator {
 public:
  std::string toString() const override { return "OutputOperator"; }
  Operator* clone() const override { return new OutputOperator(); }

  void checkAndOutput(const std::vector<CompressedSubgraphs>& input) const { /* TODO(tatiana) */
  }

  uint64_t checkAndCount(const std::vector<CompressedSubgraphs>& input) const {
    uint64_t count = 0, rough = 0;
    for (auto& group : input) {
      count += group.getNumIsomorphicSubgraphs();
      rough += group.getNumSubgraphs();
    }
    // LOG(INFO) << "# isomorphic matches " << count << " roughly " << rough;
    return count;
  }
};

}  // namespace circinus

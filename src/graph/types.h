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

#include <cinttypes>
#include <iostream>
#include <memory>
#include <vector>

namespace circinus {

using EdgeID = uint64_t;
using LabelID = uint32_t;
using QueryVertexID = uint32_t;
using VertexID = uint64_t;

using VertexSet = std::shared_ptr<std::vector<VertexID>>;

const LabelID ALL_LABEL = ~0u;

enum class GraphType : uint32_t { Normal, GraphView, BipartiteGraphView };

enum class CandidateScopeType : uint8_t { All, Partition, Inverse };

class CandidateScope {
 private:
  CandidateScopeType type_ = CandidateScopeType::All;
  uint32_t partition_ = 0;

 public:
  explicit CandidateScope(uint32_t partition) {
    partition_ = partition;
    type_ = CandidateScopeType::Partition;
  }

  void usePartition(uint32_t partition) {
    partition_ = partition;
    type_ = CandidateScopeType::Partition;
  }

  void excludePartition(uint32_t partition) {
    partition_ = partition;
    type_ = CandidateScopeType::Inverse;
  }

  inline auto getType() const { return type_; }
  inline uint32_t getPartition() const { return partition_; }

  void print(std::ostream& oss) const {
    if (type_ == CandidateScopeType::All) {
      oss << "all";
    } else {
      if (type_ == CandidateScopeType::Inverse) {
        oss << '-';
      }
      oss << partition_;
    }
  }
};

}  // namespace circinus

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
#include <cstdint>

#include "graph/types.h"
#include "utils/flags.h"

namespace circinus {
class ExecutionConfig {
 private:
  /* system level config */
  const uint32_t max_parallelism_;

  /* transient, for each operator */
  VertexID input_size_;
  uint32_t parallelism_ = 1;

 public:
  explicit ExecutionConfig(uint32_t max_parallelism = 1) : max_parallelism_(max_parallelism) {}

  uint32_t getMaxParallelism() const { return max_parallelism_; }
  uint32_t getParallelism() const { return parallelism_; }
  VertexID getInputSize() const { return input_size_; }
  uint32_t getBatchSize() const { return FLAGS_batch_size; }

  void setParallelism(uint32_t p) { parallelism_ = std::min(p, max_parallelism_); }
  void setInputSize(VertexID size) { input_size_ = size; }
};

}  // namespace circinus

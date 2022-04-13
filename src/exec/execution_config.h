#pragma once

#include <algorithm>
#include <cstdint>

#include "graph/types.h"
#include "utils/flags.h"

namespace circinus {
class ExecutionConfig {
 private:
  /* system level config */
  const uint32_t n_executors_;
  uint32_t max_parallelism_;

  /* transient, for each operator */
  VertexID input_size_;
  uint32_t parallelism_ = 1;

 public:
  explicit ExecutionConfig(uint32_t n_executors = 1, uint32_t max_parallelism = 1)
      : n_executors_(n_executors), max_parallelism_(max_parallelism) {}

  uint32_t getMaxParallelism() const { return max_parallelism_; }
  uint32_t getParallelism() const { return parallelism_; }
  VertexID getInputSize() const { return input_size_; }
  uint32_t getBatchSize() const { return FLAGS_batch_size; }
  uint32_t getNumExecutors() const { return n_executors_; }

  void setMaxParallelism(uint32_t p) { max_parallelism_ = p; }
  void setParallelism(uint32_t p) { parallelism_ = std::min(p, max_parallelism_); }
  void setInputSize(VertexID size) { input_size_ = size; }
};

}  // namespace circinus

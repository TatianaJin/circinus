#pragma once

#include <chrono>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "graph/compressed_subgraphs.h"
#include "ops/operator.h"
#include "plan/vertex_relationship.h"
#include "utils/query_utils.h"
#include "utils/utils.h"

namespace circinus {

enum class OutputType : uint32_t { Subgraph = 0, Count };

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

  inline void limit(uint64_t total_limit) { limit_per_thread_ = total_limit; }

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

  std::vector<std::pair<uint32_t, uint32_t>> set_less_than_constraints_;  // index first < index second

  /* for each same label group */
  // enumerate order of set indices when enforcing partial order
  std::vector<std::vector<uint32_t>> enumerate_orders_;
  // j: {i: the orders of smaller set indices wrt. enumerate_orders_[j][i]}
  std::vector<std::vector<std::vector<uint32_t>>> constraints_adjs_;
  const VertexRelationship* qv_relationship_ = nullptr;

 public:
  explicit OutputOperator(SameLabelIndices&& same_label_indices) : same_label_indices_(std::move(same_label_indices)) {}
  virtual ~OutputOperator() {}

  /**
   * @param same_label_indices The indices of vertices of the same label {{keys}, {sets}}
   */
  static OutputOperator* newOutputOperator(
      OutputType type, std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>&& same_label_indices);

  inline void setOutput(Outputs* outputs) { outputs_ = outputs; }
  inline Outputs* getOutput() const { return outputs_; }

  inline const auto& getSameLabelIndices() const { return same_label_indices_; }

  inline bool validateAndOutput(const std::vector<CompressedSubgraphs>& input, uint32_t output_index) const {
    uint32_t start = 0;
    return validateAndOutput(input, start, input.size(), output_index);
  }

  void setPartialOrder(std::vector<std::pair<uint32_t, uint32_t>>&& constraints);

  inline void setEquivalentClasses(const VertexRelationship& vr) { qv_relationship_ = &vr; }

  virtual bool validateAndOutput(const std::vector<CompressedSubgraphs>& input, uint32_t& input_start,
                                 uint32_t input_end, uint32_t output_index) const = 0;

  inline bool validateAndOutputAndProfile(const std::vector<CompressedSubgraphs>& input, uint32_t input_start,
                                          uint32_t input_end, uint32_t output_index, ProfileInfo* info) const {
    auto start = std::chrono::high_resolution_clock::now();
    uint32_t old_start = input_start;
    auto ret = validateAndOutput(input, input_start, input_end, output_index);
    auto stop = std::chrono::high_resolution_clock::now();
    info->total_time_in_milliseconds += toMilliseconds(start, stop);
    info->total_num_input_subgraphs += getNumSubgraphs(input, old_start, input_start);
    info->total_input_size += input_start - old_start;
    return ret;
  }

  std::string toProfileString(const ProfileInfo& info) const override {
    std::stringstream ss;
    ss << toString() << ',' << info.total_time_in_milliseconds << ',' << info.total_input_size << ",,"
       << info.total_num_input_subgraphs;
    return ss.str();
  }
};

}  // namespace circinus

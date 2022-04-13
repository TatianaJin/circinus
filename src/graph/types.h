#pragma once

#include <cinttypes>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "utils/hashmap.h"

namespace circinus {

using EdgeID = uint64_t;
using LabelID = uint32_t;
using QueryVertexID = uint32_t;
using VertexID = uint64_t;

constexpr LabelID ALL_LABEL = ~0u;
constexpr QueryVertexID DUMMY_QUERY_VERTEX = UINT32_MAX;

enum class GraphType : uint32_t { Normal, Partitioned, GraphView, BipartiteGraphView };

using PartialOrderConstraintMap =
    unordered_map<QueryVertexID, std::pair<std::vector<QueryVertexID>, std::vector<QueryVertexID>>>;

enum class CandidateScopeType : uint8_t { All = 0, Partition = 1, Range = 2, PartitionRange = 3, Inverse = 4 };
inline std::string CANDIDATE_SCOPE_TYPE_NAMES[5] = {"All", "Partition", "Range", "PRange", "Inverse"};
inline std::ostream& operator<<(std::ostream& os, CandidateScopeType type) {
  return os << CANDIDATE_SCOPE_TYPE_NAMES[(uint8_t)type];
}

class CandidateScope {
 private:
  CandidateScopeType type_ = CandidateScopeType::All;
  uint32_t partition_ = 0;
  // for range type
  uint32_t start_ = 0, end_ = 0;

 public:
  void addRange(uint32_t start, uint32_t end) {
    start_ = start;
    end_ = end;
    type_ = (type_ == CandidateScopeType::All) ? CandidateScopeType::Range : CandidateScopeType::PartitionRange;
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
  inline uint32_t getRangeStart() const { return start_; }
  inline uint32_t getRangeEnd() const { return end_; }

  void print(std::ostream& oss) const {
    oss << CANDIDATE_SCOPE_TYPE_NAMES[(uint8_t)type_];
    if ((uint8_t)type_ >> 2 & 1) {  // inverse
      oss << ' ' << '-';
    } else {
      oss << ' ';
    }
    if ((uint8_t)type_ & 1) {  // partition
      oss << partition_ << ' ';
    }
    if ((uint8_t)type_ >> 1 & 1) {  // range
      oss << '[' << start_ << ',' << end_ << "] ";
    }
  }
};

}  // namespace circinus

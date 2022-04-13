#pragma once

#include <string>
#include <utility>
#include <vector>
#include "glog/logging.h"

#include "graph/types.h"

namespace circinus {

class TreeNode {
 public:
  QueryVertexID parent_;
  uint32_t level_;
  std::vector<QueryVertexID> bn_;
  std::vector<QueryVertexID> fn_;
  std::vector<QueryVertexID> under_level_;
  std::vector<QueryVertexID> children_;
};

}  // namespace circinus

#pragma once

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/vertex_equivalence.h"
#include "graph/types.h"

namespace circinus {

class VertexRelationship {
  std::vector<QueryVertexID> set_qvs_;
  const VertexEquivalence* qv_equivalence_ = nullptr;

 public:
  explicit VertexRelationship(const VertexEquivalence& eq) : qv_equivalence_(&eq) {}

  inline bool isEquivalentByIndices(uint32_t set_index1, uint32_t set_index2) const {
    DCHECK_LT(set_index1, set_qvs_.size());
    DCHECK_LT(set_index2, set_qvs_.size());
    return qv_equivalence_->isEquivalent(set_qvs_[set_index1], set_qvs_[set_index2]);
  }

  const VertexEquivalence& getVertexEquivalence() const { return *qv_equivalence_; }

  inline void initSetQueryVertexIndices(const unordered_map<QueryVertexID, uint32_t>& set_qv_index) {
    set_qvs_.resize(set_qv_index.size());
    for (auto& p : set_qv_index) {
      set_qvs_[p.second] = p.first;
    }
  }

  std::pair<QueryVertexID, std::vector<QueryVertexID>> findReusableSet(
      QueryVertexID target, std::vector<QueryVertexID>& set_vertices,
      const unordered_set<QueryVertexID>& existing_vertices, uint64_t exclude_mask = 0) const;

  std::pair<bool, uint32_t> canBeReusedBy(QueryVertexID u1, QueryVertexID u2) const;
};

}  // namespace circinus

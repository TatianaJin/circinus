// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
// the License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

#pragma once

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

  inline void initSetQueryVertexIndices(const unordered_map<QueryVertexID, uint32_t>& set_qv_index) {
    set_qvs_.resize(set_qv_index.size());
    for (auto& p : set_qv_index) {
      set_qvs_[p.second] = p.first;
    }
  }

  std::pair<QueryVertexID, std::vector<QueryVertexID>> findReusableSet(
      QueryVertexID target, std::vector<QueryVertexID>& set_vertices,
      const unordered_set<QueryVertexID>& existing_vertices);
};

}  // namespace circinus

#pragma once

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "utils/hashmap.h"

namespace circinus {

struct PartialOrder {
  using Constraints =
      unordered_map<QueryVertexID, std::pair<unordered_set<QueryVertexID>, unordered_set<QueryVertexID>>>;

  struct EnforcePlan {
    bool has_constraint = false;
    std::vector<QueryVertexID> order;                    // i: the i-th vertex to check
    std::vector<std::vector<uint32_t>> constraints_adj;  // i: the orders of smaller vertices than `order[i]`
  };

  // i: the succinct set of query vertices smaller than query vertex i
  std::vector<std::vector<QueryVertexID>> po_constraint_adj;
  std::vector<uint32_t> order_index;  // the order of query vertex i, UINT32_MAX for qv without constraint
  std::vector<uint32_t> n_all_related_constraints;
  Constraints constraints;

  explicit PartialOrder(uint32_t n_qvs)
      : po_constraint_adj(n_qvs), order_index(n_qvs, UINT32_MAX), n_all_related_constraints(n_qvs, 0) {}

  std::ostream& printMinimum(std::ostream& stream) const;

  std::ostream& printFull(std::ostream& stream) const;

  /** Order vertices by considering only non-inducible (succinct) constraints */
  EnforcePlan orderVertices(const std::vector<QueryVertexID>& vertices) const;

  bool checkEnforcePlan(const EnforcePlan&) const;

  /** Considers all constraints related with the given vertex and return a non-inducible set of them.
   * @returns A vector of constraints represented as {less_than flag, vertex index}. */
  std::vector<std::pair<bool, uint32_t>> getConstraintsForVertex(
      QueryVertexID v, const unordered_map<QueryVertexID, uint32_t>& seen_vertices) const;

  std::vector<std::pair<bool, QueryVertexID>> getConstraintVerticesForVertex(
      QueryVertexID v, const unordered_map<QueryVertexID, uint32_t>& seen_vertices) const;

 private:
  template <bool return_qv>
  std::vector<std::pair<bool, std::conditional_t<return_qv, QueryVertexID, uint32_t>>> getConstraintsForVertexImpl(
      QueryVertexID v, const unordered_map<QueryVertexID, uint32_t>& seen_vertices) const;
};

}  // namespace circinus

#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "algorithms/partial_order.h"
#include "algorithms/vertex_equivalence.h"
#include "exec/execution_config.h"
#include "graph/graph_metadata.h"
#include "graph/types.h"
#include "ops/input_operator.h"
#include "ops/logical/compressed_input.h"
#include "plan/execution_plan.h"

namespace circinus {

class BacktrackingPlan {
  // logical
  std::vector<ExecutionPlan*> plans_;
  std::vector<std::unique_ptr<LogicalCompressedInputOperator>> input_operators_;  // size = plans_.size();
  // partitioned and parallel
  std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> partitioned_plans_;
  std::vector<bool> segment_mask_;
  const PartialOrder* partial_order_ = nullptr;
  const VertexEquivalence* qv_equivalent_classes_ = nullptr;

 public:
  std::vector<Operator*>& getOperators(uint32_t plan_idx = 0) { return plans_[plan_idx]->getOperators(); }

  inline Operator* getOutputOperator(uint32_t plan_idx = 0) { return plans_[plan_idx]->getOperators().back(); }

  std::unique_ptr<InputOperator> getInputOperator(uint32_t plan_idx = 0) {
    return input_operators_[plan_idx]->toPhysicalOperators();
  }

  inline const auto& getPlans() const { return plans_; }
  inline const auto& getPlan(uint32_t idx) const { return plans_[idx]; }
  inline auto getNumInputOperators() const { return input_operators_.size(); }

  inline uint32_t getNumPartitionedPlans() const { return partitioned_plans_.size(); }
  inline const auto& getPartitionedPlan(uint32_t idx) const { return partitioned_plans_[idx]; }

  inline bool toSegment(uint32_t idx) const {
    DCHECK_LT(idx, segment_mask_.size());
    return segment_mask_[idx];
  }

  inline void setQueryPartialOrder(const PartialOrder* po) {
    DCHECK(po != nullptr);
    partial_order_ = po;
    for (auto plan : plans_) setPartialOrderForPlan(plan);
  }

  inline void setEquivalentClasses(const VertexEquivalence& equivalent_classes) {
    qv_equivalent_classes_ = &equivalent_classes;
  }

  inline uint32_t addPlan(ExecutionPlan* plan) {
    plans_.push_back(plan);
    if (partial_order_ != nullptr) {
      setPartialOrderForPlan(plan);
    }
    return plans_.size() - 1;
  }

  /**
   * @param partitioned_plans A vector of {plan_idx, candiate scope for each query vertex}.
   */
  inline void addPartitionedPlans(std::vector<std::pair<uint32_t, std::vector<CandidateScope>>>&& partitioned_plans) {
    partitioned_plans_ = std::move(partitioned_plans);
  }

  inline void setPartitionedPlanSegmentMask(std::vector<bool>&& segment_mask) {
    segment_mask_ = std::move(segment_mask);
  }

  inline void addInputOperator(std::unique_ptr<LogicalCompressedInputOperator>&& op) {
    input_operators_.push_back(std::move(op));
  }

  inline void replaceInputOperator(uint32_t idx, std::unique_ptr<LogicalCompressedInputOperator>&& op) {
    input_operators_[idx] = std::move(op);
  }

  std::string toString() const {
    std::stringstream ss;
    for (uint32_t i = 0; i < plans_.size(); ++i) {
      ss << "[ Plan " << i << " ]\n";
      ss << "0," << input_operators_[i]->toPhysicalOperators()->toString() << std::endl;
      plans_[i]->printPhysicalPlan(ss);
    }
    if (!partitioned_plans_.empty()) {
      ss << "[ Partitions ]\n";
    }
    for (auto& plan : partitioned_plans_) {
      ss << "Plan " << plan.first;
      for (auto& scope : plan.second) {
        ss << ' ';
        scope.print(ss);
      }
      ss << '\n';
    }
    return ss.str();
  }

 private:
  void setPartialOrderForPlan(ExecutionPlan* plan) const {
    auto& po = *partial_order_;
    unordered_map<QueryVertexID, uint32_t> seen_vertices;
    seen_vertices.insert({plan->getMatchingOrder().front(), 0});
    auto op_size = plan->getOperators().size() - 1;
    // set traverse operators
    for (uint32_t i = 0; i < op_size; ++i) {
      auto traverse_op = dynamic_cast<TraverseOperator*>(plan->getOperators()[i]);
      traverse_op->setPartialOrder(po, seen_vertices);
      if (traverse_op->extend_vertex()) {
        seen_vertices.insert({traverse_op->getTargetQueryVertex(), seen_vertices.size()});
      }
    }
    // set output operator, key-to-key/key-to-set constraints should be already applied during expansion
    // output operator only needs to take care of set-to-set constraints
    auto output_op = dynamic_cast<OutputOperator*>(plan->getOperators().back());

    std::vector<std::pair<uint32_t, uint32_t>> set_constraints;
    for (QueryVertexID u = 0; u < po.po_constraint_adj.size(); ++u) {
      if (plan->isInCover(u)) continue;
      auto& smaller_qvs = po.po_constraint_adj[u];
      for (auto smaller : smaller_qvs) {
        if (!plan->isInCover(smaller)) {
          set_constraints.emplace_back(plan->getQueryVertexOutputIndex(smaller), plan->getQueryVertexOutputIndex(u));
        }
      }
    }
    output_op->setPartialOrder(std::move(set_constraints));
  }
};

}  // namespace circinus

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

#include <string>

#include "zmq.hpp"
#include "zmq_addon.hpp"

#include "exec/result.h"
#include "plan/backtracking_plan.h"

namespace circinus {

class ExecutorToCostLearnerClient {
 private:
  zmq::socket_t socket_;

 public:
  ExecutorToCostLearnerClient(zmq::context_t* zmq_ctx, const std::string& cost_learner_address)
      : socket_(*zmq_ctx, ZMQ_PUSH) {
    socket_.setsockopt(ZMQ_LINGER, 0);
    socket_.connect(cost_learner_address);
    LOG(INFO) << "Connect to cost learner at " << cost_learner_address;
  }

  void sendUpdates(ProfiledExecutionResult* result, BacktrackingPlan* plan, QueryGraph* q, uint32_t n_partitions,
                   uint32_t n_labels, const std::vector<std::vector<CandidateSetView>>& candidates) {
    CHECK_NOTNULL(q);
    CHECK_NOTNULL(result);
    CHECK_NOTNULL(plan);
    zmq::multipart_t msg;
    addGenerateQueryInfo(&msg, q, n_labels, n_partitions);
    // 2. get query vertex features for each plan
    uint32_t n_qvs = q->getNumVertices();
    auto& logical_plans = plan->getPlans();
    auto n_logical_plans = logical_plans.size();
    // assume only one partitioning qv
    auto& scopes = plan->getPartitionedPlan(0).second;
    QueryVertexID pqv = DUMMY_QUERY_VERTEX;
    for (uint32_t i = 0; i < scopes.size(); ++i) {
      if (scopes[i].getType() != CandidateScopeType::All) {
        CHECK(scopes[i].getType() == CandidateScopeType::Partition) << "Assume at most one partitioning qv";
        CHECK_EQ(pqv, DUMMY_QUERY_VERTEX) << "Assume at most one partitioning qv";
        pqv = i;
      }
    }
    // send pqv
    msg.addtyp<QueryVertexID>(pqv == DUMMY_QUERY_VERTEX ? n_qvs : pqv);
    auto n_partitioned_plans = plan->getNumPartitionedPlans();
    CHECK_EQ(candidates.size(), n_partitioned_plans);
    // 2.1 qv candidate cardinality for each logical plan
    std::vector<std::vector<uint32_t>> plan_vertex_card(n_logical_plans, std::vector<uint32_t>(n_qvs));
    for (uint32_t i = 0; i < n_partitioned_plans; ++i) {
      auto& partitioned_plan = plan->getPartitionedPlan(i);
      auto& card = plan_vertex_card[partitioned_plan.first];
      for (uint32_t v = 0; v < n_qvs; ++v) {
        if (v == pqv) {
          card[v] += candidates[i][v].size();
        } else {
          card[v] = candidates[i][v].size();
          CHECK_EQ(card[v], candidates[0][v].size());
        }
      }
    }
    // 2.2 aggregate partition scope for each logical plan
    std::vector<std::string> plan_vertex_partitions(n_logical_plans);
    for (uint32_t i = 0; i < n_partitioned_plans; ++i) {
      auto& partitioned_plan = plan->getPartitionedPlan(i);
      plan_vertex_partitions[partitioned_plan.first] = plan_vertex_partitions[partitioned_plan.first] +
                                                       std::to_string(partitioned_plan.second[pqv].getPartition()) +
                                                       ",";
    }

    msg.addtyp<uint32_t>(n_logical_plans);
    for (uint32_t i = 0; i < n_logical_plans; ++i) {
      // send card
      msg.addmem(plan_vertex_card[i].data(), n_qvs * sizeof(uint32_t));
      // send pqv partitions
      msg.addmem(plan_vertex_partitions[i].data(), plan_vertex_partitions[i].size());

      // send matching order
      auto& matching_order = logical_plans[i]->getMatchingOrder();
      msg.addmem(matching_order.data(), n_qvs * sizeof(QueryVertexID));
      auto& ops = logical_plans[i]->getOperators();
      auto& records = logical_plans[i]->getOpInputSubqueryCover();
      CHECK_EQ(records.size(), ops.size() - 1);
      msg.addtyp<uint32_t>(records.size());
      for (uint32_t opidx = 0; opidx < records.size(); ++opidx) {
        auto & [ subquery_size, cover, remove_parent_edge ] = records[opidx];
        msg.addtyp<uint32_t>(subquery_size);
        msg.addtyp<uint64_t>(cover);
        msg.addtyp<bool>(remove_parent_edge);
      }
      for (uint32_t opidx = 0; opidx < records.size(); ++opidx) {
        uint64_t truth = result->getProfile(i)[opidx + 1].total_input_size;
        msg.addtyp<uint64_t>(truth);
      }
    }
    msg.send(socket_);
  }

 private:
  inline void addGenerateQueryInfo(zmq::multipart_t* msg, QueryGraph* q, uint32_t n_labels, uint32_t n_partitions) {
    // 1. send general query graph info
    uint32_t n_qvs = q->getNumVertices();
    msg->addtyp<uint32_t>(n_labels);
    msg->addtyp<uint32_t>(n_partitions);
    msg->addtyp<uint32_t>(n_qvs);
    msg->addmem(q->getLabelArray().data(), n_qvs * sizeof(LabelID));
    for (uint32_t v = 0; v < n_qvs; ++v) {
      auto nbrs = q->getOutNeighbors(v);
      msg->addtyp<uint32_t>(nbrs.second);
      msg->addmem(nbrs.first, nbrs.second * sizeof(QueryVertexID));
    }
  }

};  // class ExecutorToCostLearnerClient

}  // namespace circinus

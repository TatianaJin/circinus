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
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "exec/result.h"
#include "graph/types.h"
#include "plan/backtracking_plan.h"
#include "plan/candidate_pruning_plan.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/query_utils.h"

namespace circinus {

/** Generates logical plans for candidate pruning and subgraph matching.
 *
 * The generated plans are to be translated into physical operators for execution.
 */
class Planner {
  QueryContext* const query_context_ = nullptr;
  CandidatePruningPlan candidate_pruning_plan_;

  // TODO(tatiana): refactor for current framework
  std::vector<std::unique_ptr<NaivePlanner>> planners_;
  std::unique_ptr<BacktrackingPlan> backtracking_plan_ = nullptr;

 public:
  explicit Planner(QueryContext& query_context) : query_context_(&query_context) {}

  virtual ~Planner() {}

  /** Generates a candidate pruning plan.
   * A candidate pruning plan has three phases, which are all optional.
   * <li> Phase 1 uses a chain of local filters to generate candidate sets for a set of query vertices
   * independently.</li>
   * <li> Phase 2 expands the generated candidate sets to the candidates of other query vertices.</li>
   * <li> Phase 3 prunes generated candidate sets all together.</li>
   */
  CandidatePruningPlan* generateCandidatePruningPlan();

  CandidatePruningPlan* updateCandidatePruningPlan(const CandidateResult* result);

  std::vector<std::vector<VertexID>> estimateCardinality() const;

  /** Generates a subgraph matching backtracking plan.
   *
   * First, for each graph partition, an execution plan containing a matching order and a compression strategy is
   * generated to respect the local graph distribution in the partition. Then the subgraph matching problem is divided
   * into multiple parallel subgraph matching tasks that find mappings of the query graph on different parts of the
   * data graph. The outputs of the parallel tasks should contain no duplicate while jointly cover all valid matches.
   *
   * Each task uses the plan generated based on a partition x, when the majority of its expanded vertices on the data
   * graph is in paritition x, so that the plan gives minimal computation workload. Since the matches may span over
   * multiple graph partitions, some tasks must work on more than one partition. To decide which plan a task should use,
   * a query vertex that has high connectivity/importance in the query graph is picked as the locality indicator. That
   * is, if a match of the indicator query vertex is in partition x, most (if not all) data vertices that appear in the
   * same embedding with this match are expected to be in partition x. The indicator serves as the primary partitioning
   * query vertex so that each task only uses a partition of the candidate set for the query vertex.
   * Then the whole search space for backtracking is divided by picking a scope of candidate vertices for each of the
   * other query vertices.
   */
  BacktrackingPlan* generateExecutionPlan(const CandidateResult*, bool multithread = true);

  BacktrackingPlan* generateExecutionPlan(const std::vector<std::vector<VertexID>>*, bool multithread = true);

 private:
  inline bool toPartitionCandidates() const {
    return query_context_->graph_metadata->numPartitions() > 1 && query_context_->query_config.use_partitioned_graph;
  }

  inline std::vector<VertexID> estimateCardinalityInner(const GraphMetadata* metadata) const {
    std::vector<VertexID> ret;
    ret.resize(query_context_->query_graph.getNumVertices(), metadata->getNumVertices());
    for (uint32_t i = 0; i < query_context_->query_graph.getNumVertices(); ++i) {
      if (metadata->hasLabelFrequency()) {
        ret[i] = std::min(ret[i], metadata->getLabelFrequency(query_context_->query_graph.getVertexLabel(i)));
      }
      if (metadata->hasOutDegreeFrequency()) {
        ret[i] = std::min(ret[i],
                          metadata->getNumVerticesWithOutDegreeGE(query_context_->query_graph.getVertexOutDegree(i)));
      }
    }
    return ret;
  }

  ExecutionPlan* generateExecutionPlan(std::vector<VertexID>& cardinality, const std::vector<QueryVertexID>& use_order,
                                       bool multithread, const std::vector<CandidateSetView>* candidate_views = nullptr,
                                       QueryVertexID partitioning_qv = UINT32_MAX);

  /* start of interface for specifying partitioning strategy */
  virtual std::vector<QueryVertexID> getPartitioningQueryVertices();  // based on query vertex occurence in covers
  std::vector<QueryVertexID> getPartitioningQueryVerticesByClosenessCentrality();
  virtual std::vector<std::pair<uint32_t, std::vector<CandidateScope>>> generatePartitionedPlans(
      const std::vector<QueryVertexID>& partitioning_qv);  // based on indicator partition only
  virtual void exhaustivePartitionPlan(std::vector<std::vector<CandidateScope>>& partitioned_plans, uint32_t level,
                                       uint32_t partition_num, const std::vector<QueryVertexID>& partitioning_qv,
                                       std::vector<CandidateScope> partitioned_plan);
  virtual void newInputOperators(const QueryGraph& q, const std::vector<QueryVertexID>& partitioning_qvs);

  virtual void newInputOperators();
  /* end of interface for specifying partitioning strategy */
};

}  // namespace circinus

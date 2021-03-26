# Copyright 2021 HDL
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

include_directories(${PROJECT_SOURCE_DIR}/src ${EXTERNAL_INCLUDES})

set(circinus_algo_sources
  ${PROJECT_SOURCE_DIR}/src/algorithms/vertex_cover.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/minimum_weight_vertex_cover.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/k_core.cc
)
set(circinus_graph_sources
  ${PROJECT_SOURCE_DIR}/src/graph/graph.cc
  ${PROJECT_SOURCE_DIR}/src/graph/query_graph.cc
)
set(circinus_ops_sources
  ${PROJECT_SOURCE_DIR}/src/ops/enumerate_key_expand_to_set_operator.cc
  ${PROJECT_SOURCE_DIR}/src/ops/expand_set_to_key_vertex_operator.cc
  ${PROJECT_SOURCE_DIR}/src/ops/expand_edge_operator.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/cfl_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/nlf_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/tso_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/dpiso_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/gql_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/order/cfl_order.cc
  ${PROJECT_SOURCE_DIR}/src/ops/order/tso_order.cc
  ${PROJECT_SOURCE_DIR}/src/ops/output_operator.cc
  ${PROJECT_SOURCE_DIR}/src/ops/scans/ldf_scan.cc
  ${PROJECT_SOURCE_DIR}/src/ops/traverse_operator.cc
)
set(circinus_plan_sources
  ${PROJECT_SOURCE_DIR}/src/plan/execution_plan.cc
  ${PROJECT_SOURCE_DIR}/src/plan/naive_planner.cc
  ${PROJECT_SOURCE_DIR}/src/plan/operator_tree.cc
)
set(circinus_utils_sources
  ${PROJECT_SOURCE_DIR}/src/utils/flags.cc
)


if (BUILD_SHARED)
  set(circinus_sources
    ${circinus_algo_sources}
    ${circinus_graph_sources}
    ${circinus_ops_sources}
    ${circinus_plan_sources}
    ${circinus_utils_sources}
    PARENT_SCOPE)
else()
  add_library(circinus_algo_obj OBJECT ${circinus_algo_sources})
  add_library(circinus_graph_obj OBJECT ${circinus_graph_sources})
  add_library(circinus_ops_obj OBJECT ${circinus_ops_sources})
  add_library(circinus_plan_obj OBJECT ${circinus_plan_sources})
  add_library(circinus_utils_obj OBJECT ${circinus_utils_sources})
  set (circinus_objs
    $<TARGET_OBJECTS:circinus_algo_obj>
    $<TARGET_OBJECTS:circinus_graph_obj>
    $<TARGET_OBJECTS:circinus_ops_obj>
    $<TARGET_OBJECTS:circinus_plan_obj>
    $<TARGET_OBJECTS:circinus_utils_obj>
    PARENT_SCOPE)
endif()

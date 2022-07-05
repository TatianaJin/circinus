include_directories(${PROJECT_SOURCE_DIR}/src ${EXTERNAL_INCLUDES})

set(circinus_algo_sources
  ${PROJECT_SOURCE_DIR}/src/algorithms/automorphism_check.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/k_core.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/leapfrog_join.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/minimum_weight_vertex_cover.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/partial_order.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/vertex_cover.cc
  ${PROJECT_SOURCE_DIR}/src/algorithms/vertex_equivalence.cc
)
set(circinus_exec_sources
  ${PROJECT_SOURCE_DIR}/src/exec/candidate_pruning_plan_driver.cc
  ${PROJECT_SOURCE_DIR}/src/exec/circinus_server.cc
  ${PROJECT_SOURCE_DIR}/src/exec/executor_manager.cc
  ${PROJECT_SOURCE_DIR}/src/exec/execution_plan_driver.cc
  ${PROJECT_SOURCE_DIR}/src/exec/result.cc
  ${PROJECT_SOURCE_DIR}/src/exec/traverse_task.cc
)
set(circinus_graph_sources
  ${PROJECT_SOURCE_DIR}/src/graph/candidate_set_view.cc
  ${PROJECT_SOURCE_DIR}/src/graph/compressed_subgraphs.cc
  ${PROJECT_SOURCE_DIR}/src/graph/graph.cc
  ${PROJECT_SOURCE_DIR}/src/graph/graph_base.cc
  ${PROJECT_SOURCE_DIR}/src/graph/graph_partition.cc
  ${PROJECT_SOURCE_DIR}/src/graph/partitioned_graph.cc
  ${PROJECT_SOURCE_DIR}/src/graph/query_graph.cc
)
set(circinus_ops_sources
  ${PROJECT_SOURCE_DIR}/src/ops/expand_edge_operator.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/gql_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/nlf_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/subgraph_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/filters/target_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/input_operator.cc
  ${PROJECT_SOURCE_DIR}/src/ops/logical/compressed_input.cc
  ${PROJECT_SOURCE_DIR}/src/ops/logical/filter/cfl_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/logical/filter/daf_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/logical/filter/gql_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/logical/filter/nlf_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/logical/filter/tso_filter.cc
  ${PROJECT_SOURCE_DIR}/src/ops/logical/scan/ldf_scan.cc
  ${PROJECT_SOURCE_DIR}/src/ops/output_operator.cc
  ${PROJECT_SOURCE_DIR}/src/ops/scans.cc
  ${PROJECT_SOURCE_DIR}/src/ops/scans/ldf_scan.cc
)
set(circinus_plan_sources
  ${PROJECT_SOURCE_DIR}/src/plan/execution_plan.cc
  ${PROJECT_SOURCE_DIR}/src/plan/naive_planner.cc
  ${PROJECT_SOURCE_DIR}/src/plan/planner.cc
  ${PROJECT_SOURCE_DIR}/src/plan/vertex_relationship.cc
)
set(circinus_utils_sources
  ${PROJECT_SOURCE_DIR}/src/utils/flags.cc
  ${PROJECT_SOURCE_DIR}/src/utils/network_utils.cc
)


if (BUILD_SHARED)
  set(circinus_sources
    ${circinus_algo_sources}
    ${circinus_exec_sources}
    ${circinus_graph_sources}
    ${circinus_ops_sources}
    ${circinus_plan_sources}
    ${circinus_utils_sources}
    PARENT_SCOPE)
else()
  add_library(circinus_algo_obj OBJECT ${circinus_algo_sources})
  add_library(circinus_exec_obj OBJECT ${circinus_exec_sources})
  add_library(circinus_graph_obj OBJECT ${circinus_graph_sources})
  add_library(circinus_ops_obj OBJECT ${circinus_ops_sources})
  add_library(circinus_plan_obj OBJECT ${circinus_plan_sources})
  add_library(circinus_utils_obj OBJECT ${circinus_utils_sources})
  set (circinus_objs
    $<TARGET_OBJECTS:circinus_algo_obj>
    $<TARGET_OBJECTS:circinus_exec_obj>
    $<TARGET_OBJECTS:circinus_graph_obj>
    $<TARGET_OBJECTS:circinus_ops_obj>
    $<TARGET_OBJECTS:circinus_plan_obj>
    $<TARGET_OBJECTS:circinus_utils_obj>
    PARENT_SCOPE)
endif()

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

find_package(GTest)

if(GTEST_FOUND)
  enable_testing()
  file(GLOB_RECURSE UNIT_TEST_SOURCE ${PROJECT_SOURCE_DIR}/src/*_unittest.cc)
  include_directories(${GTEST_INCLUDE_DIR})

  add_executable(${PROJECT_NAME}_unit_tests ${UNIT_TEST_SOURCE} unittest_main.cc)
  target_link_libraries(${PROJECT_NAME}_unit_tests circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})
  add_test(NAME AllUnitTests
    COMMAND ${PROJECT_NAME}_unit_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable(TestVertexCover unittest_main.cc test_vertex_cover_algorithm.cc test_minimum_weight_vertex_cover_algorithm.cc)
  target_link_libraries(TestVertexCover circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})
  add_test(NAME TestVertexCover
    COMMAND TestVertexCover
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable(TestFilterCandidates unittest_main.cc test_filter_candidates.cc)
  target_link_libraries(TestFilterCandidates circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})
  add_test(NAME TestFilterCandidates
    COMMAND TestFilterCandidates
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable(TestFilterAndMWVC unittest_main.cc test_filter_and_mwvc.cc)
  target_link_libraries(TestFilterAndMWVC circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})

  add_executable(TestExpandEdgeCosts unittest_main.cc test_expand_edge_costs.cc)
  target_link_libraries(TestExpandEdgeCosts circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})

  add_executable(PrintExecutionPlan unittest_main.cc print_execution_plan.cc)
  target_link_libraries(PrintExecutionPlan circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})

  add_executable(Benchmark benchmark.cc)
  target_link_libraries(Benchmark circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})

  add_executable(BenchmarkNoFilter benchmark.cc)
  target_link_libraries(BenchmarkNoFilter circinus_no_filter ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})

  add_executable(Benchmark1 benchmark1.cc)
  target_link_libraries(Benchmark1 circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})

  add_executable(DynamicCoverPreliminary dynamic_cover_preliminary.cc)
  target_link_libraries(DynamicCoverPreliminary circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY})

  add_executable(TestMetis test_metis.cc)
  target_link_libraries(TestMetis circinus ${GTEST_LIBRARY} ${EXTERNAL_LIBRARY} ${PROJECT_SOURCE_DIR}/third_party/lib/libmetis.so ${PROJECT_SOURCE_DIR}/third_party/lib/libGKlib.a )

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  add_custom_target(tests COMMAND "${CMAKE_CTEST_COMMAND}" COMMENT "Running tests..." DEPENDS ${PROJECT_NAME}_unit_tests)
else()
  message(WARNING "Google Test not found")
endif()

add_executable(TestPrintStyle test_print_style.cc)
target_link_libraries(TestPrintStyle circinus ${EXTERNAL_LIBRARY})

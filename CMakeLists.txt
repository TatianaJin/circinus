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

cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.13.0)

project(circinus VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### Options
option(BUILD_SHARED "Build shared library" ON)
option(ADDRESS_SANITIZER "Enable address sanitizer in debug mode" ON)
option(BUILD_TESTS "Build tests" ON)

# set release as the default build type
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Release, Debug, RelWithDebInfo, MinSizeRel")

# Build configuration
if(CMAKE_COMPILER_IS_GNUCXX)
  if(ADDRESS_SANITIZER) # Check memory problems
    message(STATUS "Address sanitizer is enabled in debug mode")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
  else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
  endif()
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -march=native -DNDEBUG")
endif()

### Dependencies and CMake utils
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(dep)


### Source
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/third_party/include)
add_subdirectory(src)

### Libraries
# Shared
if (BUILD_SHARED)
  add_library(circinus SHARED ${circinus_sources})
  target_compile_definitions(circinus PUBLIC "USE_FILTER")
  target_link_libraries(circinus ${EXTERNAL_LIBRARY})

  add_library(circinus_no_filter SHARED ${circinus_sources})
  target_link_libraries(circinus_no_filter ${EXTERNAL_LIBRARY})
else()
  # Static
  add_library(circinus STATIC ${circinus_objs})
  target_compile_definitions(circinus PUBLIC "USE_FILTER")
  target_link_libraries(circinus ${EXTERNAL_LIBRARY})

  add_library(circinus_no_filter STATIC ${circinus_objs})
  target_link_libraries(circinus_no_filter ${EXTERNAL_LIBRARY})
endif()
add_dependencies(circinus ${EXTERNAL_DEPENDENCY})
add_dependencies(circinus_no_filter ${EXTERNAL_DEPENDENCY})


### Executables
add_executable(Circinus src/circinus_main.cc)
target_link_libraries(Circinus circinus ${EXTERNAL_LIBRARY})

add_executable(CircinusIntersectionCache src/circinus_main.cc ${circinus_sources})
target_compile_definitions(CircinusIntersectionCache PUBLIC "USE_FILTER")
target_compile_definitions(CircinusIntersectionCache PUBLIC "INTERSECTION_CACHE")
target_link_libraries(CircinusIntersectionCache ${EXTERNAL_LIBRARY})
add_dependencies(CircinusIntersectionCache ${EXTERNAL_DEPENDENCY})

add_executable(GraphToBinary src/tools/graph_to_binary_main.cc)
target_link_libraries(GraphToBinary circinus ${EXTERNAL_LIBRARY})

# Tests
if (BUILD_TESTS)
  include(CTest)
  add_subdirectory(tests)
endif()

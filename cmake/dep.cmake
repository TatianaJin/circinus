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

set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/third_party)

### pthread
find_package(Threads)
list(APPEND EXTERNAL_LIBRARY ${CMAKE_THREAD_LIBS_INIT})

# test that filesystem header actually is there and works
try_compile(HAS_FS "${PROJECT_BINARY_DIR}/tmp"
  "${PROJECT_SOURCE_DIR}/tests/has_filesystem.cc"
  CMAKE_FLAGS -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON
  LINK_LIBRARIES stdc++fs)
if(HAS_FS)
  message(STATUS "Compiler has filesystem support")
  list(APPEND EXTERNAL_LIBRARY stdc++fs)
  list(APPEND EXTERNAL_DEFINITION "-DHAS_FILESYSTEM")
else()
  message(WARNING "Compiler is missing filesystem capabilities")
endif(HAS_FS)

### gflags
find_package(gflags REQUIRED)
message(STATUS "Found gflags:")
message(STATUS " (Headers)      ${gflags_INCLUDE_DIR}")
message(STATUS " (Library)      ${gflags_LIBRARIES}")
list(APPEND EXTERNAL_INCLUDES ${gflags_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBRARY ${gflags_LIBRARIES})


### glog
find_path(glog_INCLUDE_DIR NAMES glog/logging.h)
find_library(glog_LIBRARY NAMES glog)
if(glog_INCLUDE_DIR AND glog_LIBRARY)
  message(STATUS "Found glog:")
  message(STATUS "  (Headers)       ${glog_INCLUDE_DIR}")
  message(STATUS "  (Library)       ${glog_LIBRARY}")
  list(APPEND EXTERNAL_INCLUDES ${glog_INCLUDE_DIR})
  list(APPEND EXTERNAL_LIBRARY ${glog_LIBRARY})
else()
  message(ERROR "glog not found")
  if(INSTALL_VENDORED_LIBS)
    install(FILES ${glog_LIBRARY} DESTINATION "lib")
    install(DIRECTORY ${glog_INCLUDE_DIR}/glog DESTINATION "include")
  endif()
endif()

### gperf
set(gperf_DEFINITION "-DWITH_GPERF")
find_path(profiler_INCLUDE NAMES gperftools/profiler.h)
find_library(profiler_LIB NAMES profiler)
if(profiler_INCLUDE AND profiler_LIB)
  set(gperf_FOUND true)
  message(STATUS "Found gperf:")
  message(STATUS "  (Headers)       ${profiler_INCLUDE}")
  message(STATUS "  (Library)       ${profiler_LIB}")
  list(APPEND EXTERNAL_INCLUDES ${profiler_INCLUDE})
  list(APPEND EXTERNAL_LIBRARY ${profiler_LIB})
  list(APPEND EXTERNAL_DEFINITION ${gperf_DEFINITION})
else()
  message(STATUS "gperf not found")
endif()


add_definitions(${EXTERNAL_DEFINITION})
include_directories(${EXTERNAL_INCLUDES})

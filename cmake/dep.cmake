set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/third_party)
include(ExternalProject)

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

find_path(bliss_INCLUDE NAMES bliss/graph.hh HINTS ${THIRDPARTY_DIR})
find_library(bliss_LIBRARY NAMES bliss HINTS ${THIRDPARTY_DIR}/bliss)
message(STATUS "Finding bliss")
message(STATUS "  (Headers)    ${bliss_INCLUDE}")
message(STATUS "  (Library)    ${bliss_LIBRARY}")
if (bliss_INCLUDE AND bliss_LIBRARY)
  list(APPEND EXTERNAL_INCLUDES ${bliss_INCLUDE})
  list(APPEND EXTERNAL_LIBRARY ${bliss_LIBRARY})
else()
  message(ERROR "bliss not found")
endif()

set(install_metis false)
# Workaround for avoiding repeated build due to BUILD_IN_SOURCE
find_path(metis_INCLUDE NAMES metis.h HINTS "${THIRDPARTY_DIR}/include")
find_library(metis_LIBRARY NAMES metis HINTS ${THIRDPARTY_DIR}/lib)
if (metis_INCLUDE AND metis_LIBRARY)
  message(STATUS "Found metis")
  message(STATUS "  (Headers)     ${metis_INCLUDE}")
  message(STATUS "  (Library)     ${metis_LIBRARY}")
else()
  # Dependency for metis
  set(install_gklib false)
  find_library(gklib_LIBRARY NAMES GKlib HINTS ${THIRDPARTY_DIR}/lib)
  if (gklib_LIBRARY)
    message(STATUS "Use GKlib in ${gklib_LIBRARY}")
  else()
    set(install_gklib true)
    ###gklib
    ExternalProject_Add(gklib_ep
      CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}"
      GIT_REPOSITORY https://github.com/KarypisLab/GKlib.git
      GIT_TAG METIS-v5.1.1-DistDGL-0.5
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND make config prefix=${THIRDPARTY_DIR}
      BUILD_COMMAND make
      INSTALL_COMMAND make install
    )
  endif()

  set(install_metis true)
  ExternalProject_Add(metis_ep
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}/metis"
    GIT_REPOSITORY  https://github.com/b3ng1998/METIS.git
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND make config shared=1 i64=1 prefix=${THIRDPARTY_DIR} gklib_path=${THIRDPARTY_DIR}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  list(APPEND EXTERNAL_DEPENDENCY metis_ep)
  if (install_metis AND install_gklib)
    add_dependencies(metis_ep gklib_ep)
  endif()
endif()
list(APPEND EXTERNAL_LIBRARY ${PROJECT_SOURCE_DIR}/third_party/lib/libmetis.so)
list(APPEND EXTERNAL_LIBRARY ${PROJECT_SOURCE_DIR}/third_party/lib/libGKlib.a)


### parallel hashmap
ExternalProject_Add(phmap_ep
  CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${THIRDPARTY_DIR}"
  GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
  GIT_TAG 1.32
)
list(APPEND EXTERNAL_DEPENDENCY phmap_ep)

### zmq
find_path(ZMQ_INCLUDE_DIR NAMES zmq.hpp zmq_addon.hpp)
find_library(ZMQ_LIBRARY NAMES zmq)
if(ZMQ_INCLUDE_DIR AND ZMQ_LIBRARY)
  message(STATUS "Found ZeroMQ:")
  message(STATUS "  (Headers)       ${ZMQ_INCLUDE_DIR}")
  message(STATUS "  (Library)       ${ZMQ_LIBRARY}")
  list(APPEND EXTERNAL_INCLUDES ${ZMQ_INCLUDE_DIR})
  list(APPEND EXTERNAL_LIBRARY ${ZMQ_LIBRARY})
else()
  message(FATAL_ERROR "ZeroMQ not found")
endif()


add_definitions(${EXTERNAL_DEFINITION})
include_directories(${EXTERNAL_INCLUDES})

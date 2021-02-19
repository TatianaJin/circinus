Circinus: an efficient parallel subgraph matching framework that reduces computation redundancy and optimizes memory usage
=======

[![status](https://github.com/TatianaJin/circinus/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/TatianaJin/circinus/actions/workflows/ci.yml)

### Dependency

Circinus has the following required dependencies:

- CMake (>= 3.13.0)
- GCC (>= 8.2.0, requires c++17 support)
- Gflags
- Glog

For testing, Circinus uses [googletest](https://github.com/google/googletest/releases/tag/release-1.8.0) (needed when `cmake -DBUILD_TESTS=ON`).


### Build

1. Build and install
        $ git clone https://github.com/tatianajin/circinus.git && cd circinus
        $ mkdir -p release && cd release
        $ cmake .. -DCMAKE_BUILD_TYPE=Release # CMAKE_BUILD_TYPE: Release, Debug, RelWithDebInfo
        $ make -j4 # 4 if to use 4 threads to compile


2. Run unit tests

        $ make test                            # Run unit tests

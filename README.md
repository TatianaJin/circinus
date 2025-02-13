Circinus: an efficient parallel subgraph matching framework that reduces computation redundancy and optimizes memory usage
=======

### Disclaimer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

**This repo is no longer maintained since the owners are now working full-time, and the code is released as it is. For reproducing the results reported in the paper, please check the artifact1 branch for experiment scripts. Use at your own discretion.**

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


3. Run unit tests

        $ make test                            # Run unit tests


### Introduction to Circinus

#### Circinus Server and Client

Circinus Server manages the whole process of query processing. It owns an executor manager for query execution.

Circinus Client takes query from user, sends the query to the server, and displays the results to user. Multiple concurrent clients are supported, while a standalone mode is available for running Circius Server and Client in the same process.

#### Executor Manager

Accepts tasks from the circinus server and supports asynchronous execution of tasks. A thread pool is maintained to run tasks in parallel.

#### Planner
For each query, Circinus Server instantiates a Planner that generates logical `Plan`s for candidate pruning and subgraph matching. The generated `Plan`s are then translated into physical operators based on runtime information for execution.

A `Plan` has a `PlanDriver` that generates task intances for parallel execution of the physical operators. While the physical operators provide the application logic, the `PlanDriver` is a bridge between the execution engine and the application logic, which specifies the input to the logic, and leave task coordination to the execution engine.

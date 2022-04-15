Circinus: an efficient parallel subgraph matching framework that reduces computation redundancy and optimizes memory usage
=======

### Dependency

Circinus has the following required dependencies:

- CMake (>= 3.13.0)
- GCC (>= 8.2.0, requires c++17 support)
- Gflags
- Glog
- [bliss-0.73](http://www.tcs.hut.fi/Software/bliss/bliss-0.73.zip)

For testing, Circinus uses [googletest](https://github.com/google/googletest/releases/tag/release-1.8.0) (needed when `cmake -DBUILD_TESTS=ON`).

#### Installing Bliss

A bash script `scripts/install_bliss.sh` is provided to download and install bliss in the `third_party` folder.


### Build

1. Build and install
```bash
git clone <git link> && cd circinus
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release # CMAKE_BUILD_TYPE: Release, Debug, RelWithDebInfo
make Circinus Benckmark -j4 # 4 if to use 4 threads to compile
```

### Experimental Data and Query

The experiment results in paper is in `exp_results_in_paper/exp/data`.

If you want to reproduce the results in paper, please follow steps below.

Export the home dir of circinus.
```
export CIRCINUS_HOME=  #project dir
```

Create file to store artifact experiment results.
```
cd $CIRCINUS_HOME
mkdir -p exp_results_artifact/exp/data/
cd exp_results_artifact/exp/data/
mkdir cardinality_sensitivity  compare_peregrine  compare_single_thread  compare_unlabeled  parallelization  redundancy_multithread  vc_influence
```

All experiment scripts are in the `scripts/exp_scripts/`, check details in `scripts/exp_scripts/run_all_exps.sh`.

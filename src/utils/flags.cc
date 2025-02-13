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

#include "utils/flags.h"

#include "gflags/gflags.h"

DEFINE_int32(batch_size, 1024, "Batch size for input and output.");
DEFINE_int32(num_cores, 1, "The number of cores to use for thread pool.");
DEFINE_int32(profile, 0, "True means profiling the execution");
DEFINE_uint64(set_pruning_threshold, 0,
              "The threshold to prune by non-repeated-vertex check the same-label sets whose sizes are below the "
              "threshold. If 0, use the label frequency of the pruning query vertex. Default is 0.");

DEFINE_bool(standalone, false,
            "In standalone mode, a Circinus server is launched together with a simple command line driver.");

DEFINE_string(data_dir, "/data/share/project/haxe/data/subgraph_matching_datasets",
              "The default directory of datasets");

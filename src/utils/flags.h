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

#pragma once

#include "gflags/gflags.h"

// query execution
DECLARE_int32(batch_size);
DECLARE_int32(num_cores);
DECLARE_int32(profile);
DECLARE_uint64(set_pruning_threshold);
DECLARE_bool(label_filter);
DECLARE_int32(candidate_set_intersection);

// query optimization
DECLARE_bool(intersection_count_coefficient);
DECLARE_string(cost_learner);

DECLARE_bool(standalone);
DECLARE_string(data_dir);

DECLARE_int32(verbosity);

namespace circinus {

static constexpr int SHORT_PLANNER_LOG = 1;
static constexpr int VERBOSE_PLANNER_LOG = 3;
static constexpr int SHORT_EXECUTION_LOG = 4;
static constexpr int VERBOSE_EXECUTION_LOG = 12;
static constexpr int FULL_LOG = 15;
inline bool verbosePlannerLog() { return FLAGS_verbosity & VERBOSE_PLANNER_LOG; }
inline bool shortPlannerLog() { return FLAGS_verbosity & SHORT_PLANNER_LOG; }
inline bool shortExecutionLog() { return FLAGS_verbosity & SHORT_EXECUTION_LOG; }
inline bool verboseExecutionLog() { return FLAGS_verbosity & VERBOSE_EXECUTION_LOG; }

}  // namespace circinus

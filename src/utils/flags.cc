#include "utils/flags.h"

#include "gflags/gflags.h"

DEFINE_int32(batch_size, 1024, "Batch size for input and output.");
DEFINE_int32(num_cores, 1, "The number of cores to use for thread pool.");
DEFINE_uint64(set_pruning_threshold, 0,
              "The threshold to prune by non-repeated-vertex check the same-label sets whose sizes are below the "
              "threshold. If 0, use the label frequency of the pruning query vertex. Default is 0.");
DEFINE_bool(label_filter, true,
            "Whether to use label as a hint to filter neighbors before intersection during backtracking");
DEFINE_int32(candidate_set_intersection, 0,
             "0: infer by graph type, 1: except for last traverse op, 2: except for expand to set");

DEFINE_int32(seperate_enumeration, 0,
             "1 to split set before enumeration, 2 to seperate enumeration of sets into key in an operator. Default 0, "
             "disabled");
DEFINE_int32(task_split, 0, "0 lowest-level-first, 1 no splitting, 2 splitting on current level");

DEFINE_bool(standalone, false,
            "In standalone mode, a Circinus server is launched together with a simple command line driver.");

DEFINE_string(data_dir, "/data/share/project/haxe/data/subgraph_matching_datasets",
              "The default directory of datasets");

DEFINE_bool(intersection_count_coefficient, true,
            "If true, use estimated intersection count to calcuate costs for dynamic covers; otherwise, use estimated "
            "compressed group cardinality.");

DEFINE_string(cost_learner, "",
              "The address of the cost learner which uses machine learning to learn costs for compression plan.");

DEFINE_bool(break_symmetry, false, "If true, break symmetry to avoid automorphism in matching");

DEFINE_bool(path_card, true,
            "Whether to use path weights to estimate cardinality of cover when candidates are available");

DEFINE_int32(verbosity, circinus::FULL_LOG, "The verbosity level of logs");

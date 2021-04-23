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

#include "plan/operator_tree.h"

#include <ctime>
#include <string>
#include <vector>

#include "glog/logging.h"

#include "exec/task.h"
#include "exec/task_queue.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "ops/operators.h"
#include "utils/flags.h"
#include "utils/profiler.h"

namespace circinus {

bool OperatorTree::handleTask(Task* task, TaskQueue* queue, uint32_t thread_id) {
  auto op = operators_[task->getLevel()];

  auto traverse_op = dynamic_cast<TraverseOperator*>(op);
  if (traverse_op == nullptr) {
    auto output_op = dynamic_cast<OutputOperator*>(op);
    DCHECK(output_op != nullptr);
    return output_op->validateAndOutput(task->getInput(), thread_id);
  }

  // TODO(tatiana): data graph may change
  uint32_t last_input_index = 0;
  traverse_op->input(task->getInput(), task->getDataGraph());

  while (true) {
    std::vector<CompressedSubgraphs> outputs;
    auto start_time = clock();
    auto size = traverse_op->expand(&outputs, FLAGS_batch_size);
    if (FLAGS_profile) {
      std::string str(typeid(*traverse_op).name());
      (*profiler_)
          .addLog(task->getLevel(), str, traverse_op->getInputIndex() - last_input_index, size,
                  ((double)clock() - start_time) / CLOCKS_PER_SEC);
      last_input_index = traverse_op->getInputIndex();
    }
    if (size == 0) return false;
    queue->putTask(task->getLevel() + 1, std::move(outputs), task->getDataGraph());
  }
  return false;
}

bool OperatorTree::execute(const Graph* g, const std::vector<CompressedSubgraphs>& inputs, uint32_t level, bool useBG) {
  std::vector<CompressedSubgraphs> outputs;
  auto op = operators_[level];
  if (level == operators_.size() - 1) {
    auto output_op = dynamic_cast<OutputOperator*>(op);
    return output_op->validateAndOutput(inputs, 0);
  }
  auto traverse_op = dynamic_cast<TraverseOperator*>(op);
  uint32_t last_input_index = 0;
  traverse_op->input(inputs, g);
  if (useBG) traverse_op->useBipartiteGraph();
  while (true) {
    outputs.clear();
    auto start_time = clock();
    auto size = traverse_op->expand(&outputs, FLAGS_batch_size);
    if (size == 0) {
      break;
    }
    if (execute(g, outputs, level + 1, useBG)) {
      return true;
    }
  }
  return false;
}

bool OperatorTree::profile(const Graph* g, const std::vector<CompressedSubgraphs>& inputs, uint32_t query_type,
                           uint32_t level, bool useBG) {
  std::vector<CompressedSubgraphs> outputs;
  auto op = operators_[level];
  if (level == operators_.size() - 1) {
    auto output_op = dynamic_cast<OutputOperator*>(op);
    return output_op->validateAndOutputAndProfile(inputs, 0);
  }
  auto traverse_op = dynamic_cast<TraverseOperator*>(op);
  uint32_t last_input_index = 0;
  traverse_op->inputAndProfile(inputs, g);
  if (useBG) traverse_op->useBipartiteGraph();
  while (true) {
    outputs.clear();
    auto start_time = clock();
    auto size = traverse_op->expandAndProfile(&outputs, FLAGS_batch_size, query_type);
    if (size == 0) {
      break;
    }
    if (profile(g, outputs, query_type, level + 1, useBG)) {
      return true;
    }
  }
  return false;
}

}  // namespace circinus

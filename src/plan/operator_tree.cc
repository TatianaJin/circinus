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

#include <vector>

#include "glog/logging.h"

#include "exec/task.h"
#include "exec/task_queue.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "ops/operators.h"
#include "utils/flags.h"

namespace circinus {

bool OperatorTree::handleTask(Task* task, TaskQueue* queue, uint32_t thread_id) {
  auto op = operators_[task->getLevel()];

  // TODO(tatiana): handle non-traverse operators
  auto traverse_op = dynamic_cast<TraverseOperator*>(op);
  if (traverse_op == nullptr) {
    auto output_op = dynamic_cast<OutputOperator*>(op);
    DCHECK(output_op != nullptr);
    return output_op->validateAndOutput(task->getInput(), thread_id);
  }

  // TODO(tatiana): data graph may change
  traverse_op->input(task->getInput(), task->getDataGraph());

  while (true) {
    std::vector<CompressedSubgraphs> outputs;
    auto size = traverse_op->expand(&outputs, FLAGS_batch_size);
    if (size == 0) return false;
    queue->putTask(task->getLevel() + 1, std::move(outputs), task->getDataGraph());
  }
  return false;
}

bool OperatorTree::handleInput(const Graph* g, const std::vector<CompressedSubgraphs>& inputs, uint32_t level) {
  std::vector<CompressedSubgraphs> outputs;
  auto op = operators_[level];
  if (level == operators_.size() - 1) {
    auto output_op = dynamic_cast<OutputOperator*>(op);
    return output_op->validateAndOutput(inputs, 0);
  }
  auto traverse_op = dynamic_cast<TraverseOperator*>(op);
  traverse_op->input(inputs, g);
  while (true) {
    outputs.clear();
    auto size = traverse_op->expand(&outputs, FLAGS_batch_size);
    if (size == 0) {
      break;
    }
    if (handleInput(g, outputs, level + 1)) {
      return true;
    }
  }
  return false;
}

}  // namespace circinus

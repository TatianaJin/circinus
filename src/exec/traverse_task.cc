// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
// the License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

#include "exec/traverse_task.h"

namespace circinus {

void TraverseChainTask::run(uint32_t executor_idx) {
  if (task_status_ == TaskStatus::Normal) {
    if (input_size_ == 0) {
      if (verboseExecutionLog()) {
        std::stringstream ss;
        for (auto& x : *candidates_) {
          ss << x.size() << " ";
        }
        LOG(INFO) << "-------- Task " << task_id_ << " candidate sizes " << ss.str();
      }
      inputs_ = std::move(input_op_->getInputs(graph_, *candidates_));
      // if input is pruned to empty, return
      if (inputs_.empty()) {
        return;
      }
      input_size_ = inputs_.size();
    }
    setupOutputs();
    if (!setupTraverseContexts()) return;
  }
  uint64_t old_count = 0;
  if (end_level_ == operators_->size() - 1) {
    old_count = dynamic_cast<OutputOperator*>(operators_->back())->getOutput()->getCount(executor_idx);
  }
  // TODO(engineering): support match limit
  execute<QueryType::Execute>(inputs_, start_index_, input_size_, start_level_, executor_idx);
  checkSplits();
  logTask(old_count, executor_idx);
}

void TraverseChainTask::profile(uint32_t executor_idx) {
  if (task_status_ == TaskStatus::Normal) {
    profile_info_.resize(operators_->size() + 1);  // traverse chain + input operator
    if (input_size_ == 0) {
      input_op_->inputAndProfile(graph_, *candidates_, &inputs_, &profile_info_.front());
      if (inputs_.empty()) return;
      input_size_ = inputs_.size();
    }
    setupOutputs();
    if (!setupTraverseContexts()) return;
  }
  uint64_t old_count = 0;
  if (end_level_ == operators_->size() - 1) {
    old_count = dynamic_cast<OutputOperator*>(operators_->back())->getOutput()->getCount(executor_idx);
  }
  // TODO(engineering): support match limit
  execute<QueryType::Profile>(inputs_, start_index_, input_size_, start_level_, executor_idx);
  if (task_status_ == TaskStatus::Normal) {
    DCHECK_LT(start_level_, operators_->size());
    DCHECK_LT(end_level_, operators_->size());
    for (uint32_t i = start_level_; i < end_level_; ++i) {
      profile_info_[i + 1] += *traverse_context_[i - start_level_];
    }
    checkSplits();
  }
  logTask(old_count, executor_idx);
}

void TraverseChainTask::logTask(uint64_t old_count, uint32_t executor_idx) {
  if (shortExecutionLog()) {
    uint64_t new_count = 0;
    if (end_level_ == operators_->size() - 1) {
      new_count = dynamic_cast<OutputOperator*>(operators_->back())->getOutput()->getCount(executor_idx);
    }
    auto end = std::chrono::steady_clock::now();
    if (task_status_ == TaskStatus::Suspended || task_status_ == TaskStatus::Split) {
      LOG(INFO) << "Suspended Task " << task_id_ << " split " << split_level_ << " size " << splits_.size()
                << " suspend " << suspended_level_ << '/' << (operators_->size() - 1)
                << " time usage: " << toSeconds(start_time_, end) << "s. suspend interval "
                << (suspend_interval_ == nullptr ? 0 : *suspend_interval_);
    } else if (toSeconds(start_time_, end) > 1) {
      if (start_level_ == 0) {
        LOG(INFO) << "Task " << task_id_ << " input " << inputs_.size() << '/'
                  << getNumSubgraphs(inputs_, 0, inputs_.size()) << " count " << (new_count - old_count) << " ("
                  << old_count << " to " << new_count << "), time usage: " << toSeconds(start_time_, end)
                  << "s. suspend interval " << (suspend_interval_ == nullptr ? 0 : *suspend_interval_);
      } else {
        LOG(INFO) << "Inherit Task " << task_id_ << ':' << start_level_ << " input " << inputs_.size() << " count "
                  << (new_count - old_count) << " (" << old_count << " to " << new_count
                  << "), time usage: " << toSeconds(start_time_, end) << "s.";
      }
    }
  }
}

template <QueryType mode>
bool TraverseChainTask::splitInput(bool split_on_suspended_level) {
  if (split_level_ <= start_level_) {  // if split from start level, directly get from input
    if (traverse_context_[0]->canSplitInput()) {
      auto[ptr, size] = traverse_context_[0]->splitInput();
      // LOG(INFO) << "Split start " << start_level_ << " size " << size << " interval " << suspend_interval_;
      while (size > batch_size_) {
        splits_.emplace_back(start_level_, std::vector<CompressedSubgraphs>(ptr, ptr + batch_size_));
        ptr += batch_size_;
        size -= batch_size_;
      }
      splits_.emplace_back(start_level_, std::vector<CompressedSubgraphs>(ptr, ptr + size));
      if (splits_.size() == split_size_) {
        return true;
      }
    }
    split_level_ = start_level_ + 1;
  }
  auto split_max_level = std::min(suspended_level_ + split_on_suspended_level, (uint32_t)operators_->size() - 1);
  if (split_level_ == split_max_level) return !splits_.empty() && split_level_ >= operators_->size() - 2;

  if (traverse_context_[split_level_ - start_level_]->canSplitInput()) {
    auto[ptr, size] = traverse_context_[split_level_ - start_level_]->splitInput();
    // LOG(INFO) << "Split input " << split_level_ << " size " << size << " interval " << suspend_interval_;
    splits_.emplace_back(split_level_, std::vector<CompressedSubgraphs>(ptr, ptr + size));
    if (splits_.size() == split_size_) {
      return true;
    }
  }

  // split from intermediate level
  while (split_level_ < split_max_level && splits_.size() < split_size_) {
    auto ctx = traverse_context_[split_level_ - start_level_ - 1].get();
    auto original_buffer = ctx->getOutputs();
    std::vector<CompressedSubgraphs> temp_output_buffer;
    // auto split_batch = (1u << split_level_);
    // auto split_batch_size = split_batch * batch_size_;
    auto split_batch = split_level_;
    auto split_batch_size = split_level_ * batch_size_;
    temp_output_buffer.reserve(split_batch_size);
    for (uint32_t i = 0; i < split_batch; ++i) {
      temp_output_buffer.insert(temp_output_buffer.end(), original_buffer->begin(), original_buffer->end());
    }
    auto original_size = ctx->getOutputSize();
    ctx->setOutputBuffer(temp_output_buffer, 0);
    auto traverse_op = dynamic_cast<TraverseOperator*>((*operators_)[split_level_ - 1]);
    uint32_t size = 0;
    if
      constexpr(isProfileMode(mode)) size = traverse_op->expandAndProfile(split_batch_size, ctx);
    else
      size = traverse_op->expand(split_batch_size, ctx);
    ctx->setOutputBuffer(*original_buffer, original_size);
    if (size == 0) {
      ++split_level_;
    } else {
      temp_output_buffer.erase(temp_output_buffer.begin() + size, temp_output_buffer.end());
      // LOG(INFO) << "Split middle " << split_level_ << " size " << size << " interval " << suspend_interval_;
      splits_.emplace_back(split_level_, std::move(temp_output_buffer));
    }
  }
  return splits_.size() == split_size_ || (!splits_.empty() && split_level_ >= operators_->size() - 2);
}

template <QueryType mode>
bool TraverseChainTask::execute(const std::vector<CompressedSubgraphs>& input, uint32_t start_index,
                                uint32_t input_size, uint32_t level, uint32_t executor_idx) {
  if (isTimeOut()) {
    return true;
  }
  auto op = (*operators_)[level];
  if (level == operators_->size() - 1) {
    auto output_op = dynamic_cast<OutputOperator*>(op);
    if
      constexpr(isProfileMode(mode)) {
        return output_op->validateAndOutputAndProfile(input, 0, input_size, executor_idx, &profile_info_[level + 1]);
      }
    uint32_t start = 0;
    return output_op->validateAndOutput(input, start, input_size, executor_idx);
  }

  auto traverse_op = dynamic_cast<TraverseOperator*>(op);
  auto ctx = traverse_context_[level - start_level_].get();
  if (task_status_ == TaskStatus::Normal) {
    DCHECK_NOTNULL(ctx->getOutputs());
    ctx->setInput(input, start_index, start_index + input_size);
    if (suspend_interval_ != nullptr && *suspend_interval_ > 0 &&
        toSeconds(start_time_, std::chrono::steady_clock::now()) >= *suspend_interval_) {
      suspended_level_ = level;
      if (splitInput<mode>(traverse_op->enumeratesSet())) {  // do not suspend if no splits
        task_status_ = TaskStatus::Suspended;
        return true;
      }
    }
  } else if (level == suspended_level_) {  // when a suspended task is resumed, we need to first recreate the stack
                                           // until the last level
    task_status_ = TaskStatus::Normal;
  }

  bool finished = false;
  while (true) {
    if (isTimeOut()) {
      return true;
    }
    uint32_t size = 0;
    if (task_status_ == TaskStatus::Normal) {
      ctx->resetOutputs();
      DCHECK_EQ(ctx->getOutputs()->size(), batch_size_) << "level " << level;
      if
        constexpr(isProfileMode(mode)) size = traverse_op->expandAndProfile(batch_size_, ctx);
      else
        size = traverse_op->expand(batch_size_, ctx);

      if (size == 0) {
        break;
      }
    }
    DCHECK_EQ(outputs_[level - start_level_].size(), batch_size_) << "level " << level;
    if (execute<mode>(outputs_[level - start_level_], 0, size, level + 1, executor_idx)) {
      finished = true;
      break;
    }
  }
  if
    constexpr(isProfileMode(mode)) {
      if (task_status_ == TaskStatus::Normal) {
        ctx->total_input_size += ctx->getTotalInputSize();
      }
    }

  return finished;
}

bool TraverseChainTask::setupTraverseContexts() {
  if (candidates_ == nullptr || FLAGS_candidate_set_intersection == 3) return setupTraverseContextsWithoutCandidate();
  auto op = (*operators_)[start_level_];
  auto traverse = dynamic_cast<TraverseOperator*>(op);
  traverse_context_.reserve(end_level_ - start_level_);
  for (uint32_t level = start_level_; level < end_level_; ++level) {
    DCHECK(traverse != nullptr);
    DCHECK_LT(traverse->getTargetQueryVertex(), candidates_->size());
    const unordered_set<VertexID>* hashmap = nullptr;
    if (candidate_hashmaps_ != nullptr) {
      DCHECK_LT(traverse->getTargetQueryVertex(), candidate_hashmaps_->size());
      hashmap = (*candidate_hashmaps_)[traverse->getTargetQueryVertex()].get();
    }
    // now assume all query vertices have candidate sets
    auto ctx = createTraverseContext(&(*candidates_)[traverse->getTargetQueryVertex()], outputs_[level - start_level_],
                                     level, traverse, query_type_, hashmap);
    if (ctx == nullptr) return false;
    traverse_context_.emplace_back(std::move(ctx));
    traverse = dynamic_cast<TraverseOperator*>(traverse->getNext());
  }
  return true;
}

bool TraverseChainTask::setupTraverseContextsWithoutCandidate() {
  auto op = (*operators_)[start_level_];
  auto traverse = dynamic_cast<TraverseOperator*>(op);
  traverse_context_.reserve(end_level_ - start_level_);
  for (uint32_t level = start_level_; level < end_level_; ++level) {
    DCHECK(traverse != nullptr);
    auto ctx = createTraverseContext(nullptr, outputs_[level - start_level_], level, traverse, query_type_, nullptr);
    if (ctx == nullptr) return false;
    traverse_context_.emplace_back(std::move(ctx));
    traverse = dynamic_cast<TraverseOperator*>(traverse->getNext());
  }
  return true;
}

void TraverseChainTask::setupOutputs() {
  uint32_t key_size = inputs_.front().getKeys().size();
  uint32_t set_size = inputs_.front().getSets().size();

  auto output_size = std::make_pair(key_size, key_size + set_size);
  outputs_.resize(end_level_ - start_level_);
  for (uint32_t i = start_level_; i < end_level_; ++i) {
    CHECK_LT(i, operators_->size());
    auto op = operators_->at(i);
    output_size = ((const TraverseOperator*)op)->getOutputSize(output_size);
    outputs_[i - start_level_].resize(batch_size_, CompressedSubgraphs(output_size.first, output_size.second));
  }
}

std::vector<GraphView<GraphPartitionBase>> TraverseTask::setupGraphView(const ReorderedPartitionedGraph* g,
                                                                        const std::vector<Operator*>& operators,
                                                                        const std::vector<CandidateScope>& scopes) {
  std::vector<GraphView<GraphPartitionBase>> data_graphs_for_operators;
  CHECK_GT(operators.size(), 0);
  data_graphs_for_operators.reserve(end_level_ - start_level_);
  for (size_t i = start_level_; i < end_level_; ++i) {
    auto op = operators[i];
    auto traverse_op = dynamic_cast<TraverseOperator*>(op);
    CHECK(traverse_op != nullptr) << i << '/' << end_level_ << ' ' << op->toString();
    data_graphs_for_operators.emplace_back(traverse_op->computeGraphPartitions(g, scopes));
  }
  return data_graphs_for_operators;
}

}  // namespace circinus

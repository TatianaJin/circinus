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

#include <cinttypes>
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "graph/compressed_subgraphs.h"

namespace circinus {

struct ProfileContent {
  uint32_t level_;
  std::string op_name_;
  uint32_t input_size_;
  uint32_t output_size_;
  double time_usage_;  // s

  ProfileContent& operator+=(const ProfileContent& pc) {
    this->level_ = pc.level_;
    this->op_name_ = pc.op_name_;
    this->input_size_ += pc.input_size_;
    this->output_size_ += pc.output_size_;
    this->time_usage_ += pc.time_usage_;
    return *this;
  }
};

class Profiler {
 public:
  ~Profiler() {}

  void addLog(uint32_t level, std::string& op_name, uint32_t input_size, uint32_t output_size, double time_usage) {
    mtx_.lock();
    logs_.emplace_back(ProfileContent{level, op_name, input_size, output_size, time_usage});
    mtx_.unlock();
  }

  void addInput(const std::vector<CompressedSubgraphs> inputs, uint32_t start_idx, uint32_t end_idx) {
    mtx_.lock();
    std::vector<CompressedSubgraphs> ins;
    ins.reserve(end_idx - start_idx);
    for (uint32_t i = start_idx; i < end_idx; ++i) {
      ins.emplace_back(inputs[i]);
    }
    op_inputs_.push_back(std::move(ins));
    mtx_.unlock();
  }

  void profile(std::ostream* out) {
    (*out) << "op_name,input_size,output_size,time_usage,time_usage_/input_size_,time_usage_/output_size_\n";
    std::map<uint32_t, ProfileContent> total;
    uint32_t idx = 0;
    for (auto& pc : logs_) {
      (*out) << pc.level_ << ',' << pc.op_name_ << ',' << pc.input_size_ << ',' << pc.output_size_ << ','
             << pc.time_usage_ << ',' << pc.time_usage_ / pc.input_size_ << ',' << pc.time_usage_ / pc.output_size_
             << '\n';
      if (op_inputs_.size() != 0) {
        for (auto& in : op_inputs_[idx]) {
          (*out) << in.toString() << ",";
        }
        idx++;
        (*out) << "\n";
      }
      // group by op name
      total[pc.level_] += pc;
    }

    (*out) << "------------------------------------------------------------------------------------------\n";
    for (auto& it : total) {
      auto& pc = it.second;
      (*out) << it.first << ',' << pc.op_name_ << ',' << pc.input_size_ << ',' << pc.output_size_ << ','
             << pc.time_usage_ << ',' << pc.time_usage_ / pc.input_size_ << ',' << pc.time_usage_ / pc.output_size_
             << '\n';
    }
  }

 private:
  std::mutex mtx_;
  std::vector<ProfileContent> logs_;
  std::vector<std::vector<CompressedSubgraphs>> op_inputs_;
};

}  // namespace circinus

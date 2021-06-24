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

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "graph/candidate_set_view.h"

using circinus::VertexID;
using circinus::CandidateScope;
using circinus::CandidateSetView;

class TestCandidateSetView : public testing::Test {
 protected:
  static std::vector<VertexID> generateCandidateSet(uint32_t size) {
    std::vector<VertexID> ret(size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, size * 10);
    for (uint32_t i = 0; i < size; ++i) {
      ret[i] = distrib(gen);
    }
    std::sort(ret.begin(), ret.end());
    return ret;
  }

  static std::string toString(const std::vector<VertexID>& vec) {
    std::stringstream ss;
    ss << "size=" << vec.size();
    for (auto v : vec) {
      ss << ' ' << v;
    }
    return ss.str();
  }
};

TEST_F(TestCandidateSetView, Iterable) {
  auto candidate_set = generateCandidateSet(1000);
  std::vector<VertexID> partition_offsets = {0, 200, 500, 500, 1000};  // 4 partitions
  partition_offsets.back() = candidate_set.size();

  CandidateScope scope;
  {
    // scope: all
    CandidateSetView view(&candidate_set, scope, partition_offsets);
    uint32_t idx = 0;
    for (auto v : view) {
      ASSERT_EQ(candidate_set[idx], v) << idx << " candidate " << toString(candidate_set) << "\n partition "
                                       << toString(partition_offsets);
      ++idx;
    }
    ASSERT_EQ(idx, 1000);
  }
  {
    // scope: non-empty partition
    scope.usePartition(1);
    CandidateSetView view(&candidate_set, scope, partition_offsets);
    uint32_t idx = 200;
    for (auto v : view) {
      ASSERT_EQ(candidate_set[idx], v) << idx << " candidate " << toString(candidate_set) << "\n partition "
                                       << toString(partition_offsets);
      ++idx;
    }
    ASSERT_EQ(idx, 500);
  }
  {
    // scope: empty partition
    scope.usePartition(2);
    CandidateSetView view(&candidate_set, scope, partition_offsets);
    EXPECT_TRUE(view.empty());
    for (auto v : view) {
      ASSERT_TRUE(false) << "The view should be empty";
    }
  }
  {
    // scope: exclude non-empty partition
    scope.excludePartition(0);
    CandidateSetView view(&candidate_set, scope, partition_offsets);
    uint32_t idx = 200;
    for (auto v : view) {
      ASSERT_EQ(candidate_set[idx], v) << idx << " candidate " << toString(candidate_set) << "\n partition "
                                       << toString(partition_offsets);
      ++idx;
    }
    ASSERT_EQ(idx, 1000);
  }
  {
    // scope: exclude non-empty partition
    scope.excludePartition(0);
    CandidateSetView view(&candidate_set, scope, partition_offsets);
    uint32_t idx = 200;
    for (auto v : view) {
      ASSERT_EQ(candidate_set[idx], v) << idx << " candidate " << toString(candidate_set) << "\n partition "
                                       << toString(partition_offsets);
      ++idx;
    }
    ASSERT_EQ(idx, 1000);
  }
  {
    // scope: exclude empty partition
    scope.excludePartition(2);
    CandidateSetView view(&candidate_set, scope, partition_offsets);
    uint32_t idx = 0;
    for (auto v : view) {
      ASSERT_EQ(candidate_set[idx], v) << idx << " candidate " << toString(candidate_set) << "\n partition "
                                       << toString(partition_offsets);
      ++idx;
    }
    ASSERT_EQ(idx, 1000);
  }
}

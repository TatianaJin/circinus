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

#include <errno.h>
#include <pthread.h>
#include <thread>
#include <vector>

namespace circinus {

inline bool pinToCores(std::vector<std::thread>& threads) {
  cpu_set_t cpuset;
  for (uint32_t cpu_idx = 0; cpu_idx < threads.size(); ++cpu_idx) {
    CPU_ZERO(&cpuset);
    CPU_SET(cpu_idx, &cpuset);
    int rc = pthread_setaffinity_np(threads[cpu_idx].native_handle(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
      return false;
    }
  }
  return true;
}

}  // namespace circinus

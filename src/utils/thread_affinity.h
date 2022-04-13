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

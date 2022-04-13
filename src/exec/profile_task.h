#pragma once

#include <utility>

#include "exec/task.h"

namespace circinus {

template <typename T>
class ProfileTask : public T {
  static_assert(std::is_base_of<TaskBase, T>::value, "T should inherit from TaskBase");

 public:
  template <typename... Args>
  ProfileTask(QueryId qid, TaskId tid, std::chrono::time_point<std::chrono::steady_clock> stop_time, Args&&... args)
      : T(qid, tid, stop_time, std::forward<Args>(args)...) {}

  void run(uint32_t executor_idx) override { T::profile(executor_idx); }

  const GraphBase* getDataGraph() const override { return T::getDataGraph(); }
};

}  // namespace circinus

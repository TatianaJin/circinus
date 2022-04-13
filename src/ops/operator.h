#pragma once

#include <string>

#include "utils/utils.h"

namespace circinus {

class ProfileInfo;  // forward declaration

class Operator {
 private:
  Operator* next_ = nullptr;

 protected:
  uint32_t parallelism_;

 public:
  explicit Operator(uint32_t parallelism = 1) : parallelism_(parallelism) {}
  virtual ~Operator() {}

  inline void setNext(Operator* next) { next_ = next; }
  inline Operator* getNext() const { return next_; }

  virtual std::string toString() const { return getTypename(*this); }

  virtual std::string toProfileString(const ProfileInfo&) const { return toString(); }

  inline uint32_t getParallelism() const { return parallelism_; }

  inline void setParallelism(uint32_t parallelism) { parallelism_ = parallelism; }
};

}  // namespace circinus

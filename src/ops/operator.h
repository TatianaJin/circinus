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

#include <string>

#include "utils/utils.h"

namespace circinus {

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

  virtual std::string toProfileString() const { return toString(); }

  inline uint32_t getParallelism() const { return parallelism_; }

  inline void setParallelism(uint32_t parallelism) { parallelism_ = parallelism; }
};

}  // namespace circinus

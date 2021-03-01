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

namespace circinus {

class Operator {
 private:
  Operator* next_ = nullptr;

 public:
  virtual ~Operator() {}

  inline void setNext(Operator* next) { next_ = next; }
  inline Operator* getNext() const { return next_; }

  virtual std::string toString() const { return "Operator"; }

  // TODO(tatiana): profile info
  virtual std::string toProfileString() const { return "Operator"; }

  virtual Operator* clone() const = 0;
};

}  // namespace circinus

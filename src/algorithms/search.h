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

#define BINARY_SEARCH_THRESHOLD 32

#include <algorithm>
#include <type_traits>

namespace circinus {

template <typename Iter>
class HasLowerBound {
 private:
  template <typename U>
  static constexpr auto check(int) -> typename std::is_same<
      decltype(std::declval<U>().getLowerBound(std::declval<U>(), std::declval<typename U::value_type>())), U>::type;

  template <typename>
  static constexpr std::false_type check(...);

 public:
  using type = decltype(check<Iter>(0));
  static constexpr bool value = type::value;
};

template <typename ForwardIt, typename T>
inline constexpr ForwardIt lowerBoundLinear(ForwardIt first, ForwardIt last, const T& value) {
  while (first != last && *first < value) {
    ++first;
  }
  return first;
}

template <typename ForwardIt, typename T>
inline constexpr typename std::enable_if<HasLowerBound<ForwardIt>::value, ForwardIt>::type lowerBound(ForwardIt first,
                                                                                                      ForwardIt last,
                                                                                                      const T& value) {
  return first.getLowerBound(last, value);
}

template <typename ForwardIt, typename T>
inline constexpr typename std::enable_if<!HasLowerBound<ForwardIt>::value, ForwardIt>::type lowerBound(ForwardIt first,
                                                                                                       ForwardIt last,
                                                                                                       const T& value) {
  if (std::distance(first, last) >= BINARY_SEARCH_THRESHOLD) {
    return std::lower_bound(first, last, value);
  }
  return lowerBoundLinear(first, last, value);
}

}  // namespace circinus

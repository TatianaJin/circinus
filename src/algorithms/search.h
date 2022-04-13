#pragma once

#define BINARY_SEARCH_THRESHOLD 256

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

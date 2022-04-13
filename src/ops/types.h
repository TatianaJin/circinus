#pragma once

#include <cinttypes>

#include "utils/query_utils.h"

namespace circinus {

inline constexpr bool isProfileMode(QueryType profile) {
  return profile == QueryType::Profile || profile == QueryType::ProfileWithMiniIntersection ||
         profile == QueryType::ProfileCandidateSIEffect;
}

inline constexpr bool isProfileWithMiniIntersectionMode(QueryType profile) {
  return profile == QueryType::ProfileWithMiniIntersection;
}

inline constexpr bool isProfileCandidateSIEffect(QueryType profile) {
  return profile == QueryType::ProfileCandidateSIEffect;
}

}  // namespace circinus

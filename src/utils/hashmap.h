#pragma once

#include <map>
#include <unordered_map>
#include <unordered_set>

#ifdef USE_STL

namespace circinus {

template <typename T>
using unordered_set = std::unordered_set<T>;

template <typename K, typename V>
using unordered_map = std::unordered_map<K, V>;

template <typename K, typename V>
using ordered_map = std::map<K, V>;

}  // namespace circinus

#else  // USE_STL

#include "parallel_hashmap/btree.h"
#include "parallel_hashmap/phmap.h"

namespace circinus {

template <typename T>
using unordered_set = phmap::flat_hash_set<T>;

template <typename K, typename V>
using unordered_map = phmap::flat_hash_map<K, V>;

template <typename K, typename V>
using ordered_map = phmap::btree_map<K, V>;

}  // namespace circinus

#endif  // USE_STL

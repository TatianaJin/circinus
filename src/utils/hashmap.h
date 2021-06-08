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

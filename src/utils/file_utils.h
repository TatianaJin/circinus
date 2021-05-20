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

#ifdef HAS_FILESYSTEM
#include <filesystem>
#endif

#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

namespace circinus {

#ifdef HAS_FILESYSTEM

class Path {
 public:
  template <typename... Args>
  static inline std::string join(const std::string& dir, Args... args) {
    std::filesystem::path path(dir);
    return join(path, args...).string();
  }

  static inline bool isRelative(const std::string& path) { return std::filesystem::path(path).is_relative(); }

 private:
  static inline std::filesystem::path& join(std::filesystem::path& path, std::string& append) { return path /= append; }
  template <typename... Args>
  static inline std::filesystem::path& join(std::filesystem::path& path, std::string& append, Args... args) {
    path /= append;
    return join(path, args...);
  }
};

#else  // HAS_FILESYSTEM

/** Now only supports Linux */
class Path {
 public:
  template <typename... Args>
  static inline std::string join(const std::string& dir, Args... args) {
    std::stringstream ss;
    ss << dir;
    return join(ss, args...).str();
  }

 private:
  static inline std::stringstream& join(std::stringstream& ss, const std::string& base) { return ss << "/" << base; }
  template <typename... Args>
  static inline std::stringstream& join(std::stringstream& ss, const std::string& base, Args... args) {
    ss << "/" << base;
    return join(ss, args...);
  }
};

#endif  // HAS_FILESYSTEM

inline std::ifstream openFile(const std::string& path) {
  std::ifstream infile(path);
  if (!infile.is_open()) {
    std::stringstream ss;
    ss << std::strerror(errno) << ": " << path;
    throw std::runtime_error(ss.str());
  }
  return infile;
}

}  // namespace circinus

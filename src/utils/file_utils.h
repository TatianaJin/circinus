#pragma once

#if __GNUC__ > 7
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
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
  static inline std::string join(const std::string& dir, const std::string& base, Args... args) {
    fs::path path(dir);
    return join(path, base, args...).string();
  }

  static inline bool isRelative(const std::string& path) { return fs::path(path).is_relative(); }

 private:
  static inline fs::path& join(fs::path& path, const std::string& append) { return path /= append; }
  template <typename... Args>
  static inline fs::path& join(fs::path& path, const std::string& append, Args... args) {
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

  static inline bool isRelative(const std::string& path) { return path.empty() || path[0] != '/'; }

 private:
  static inline std::stringstream& join(std::stringstream& ss, const std::string& base) { return ss << "/" << base; }
  template <typename... Args>
  static inline std::stringstream& join(std::stringstream& ss, const std::string& base, Args... args) {
    ss << "/" << base;
    return join(ss, args...);
  }
};

#endif  // HAS_FILESYSTEM

inline std::ifstream openFile(const std::string& path, std::ios::openmode openmode = std::ios::in) {
  std::ifstream infile(path, openmode);
  if (!infile.is_open()) {
    std::stringstream ss;
    ss << std::strerror(errno) << ": " << path;
    throw std::runtime_error(ss.str());
  }
  infile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  return infile;
}

inline std::ofstream openOutputFile(const std::string& path, std::ios::openmode openmode = std::ios::out) {
  std::ofstream outfile(path, openmode);
  if (!outfile.is_open()) {
    std::stringstream ss;
    ss << std::strerror(errno) << ": " << path;
    throw std::runtime_error(ss.str());
  }
  outfile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  return outfile;
}

}  // namespace circinus

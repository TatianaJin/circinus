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
#include <unordered_map>
#include <utility>
#include <vector>

#include "graph/graph_base.h"
#include "graph/graph_metadata.h"
#include "graph/query_graph.h"
#include "graph/types.h"

#include "utils/hashmap.h"

namespace circinus {

using QueryId = uint16_t;
using TaskId = uint16_t;

enum class QueryType : uint8_t { Execute = 0, Profile, ProfileWithMiniIntersection, SampleExecute };

struct ServerEvent {
  enum Type : uint32_t {
    /* User commands */
    ShutDown = 0,
    LoadGraph,
    NewQuery,
    Explain,
    Profile,
    /* Internal phases */
    CandidatePhase,
    ExecutionPhase,
  } type;
  std::vector<std::string> args;
  std::string client_addr;
  void* data;
  QueryId query_id;

  static inline const std::vector<std::string> type_names = {"ShutDown", "LoadGraph",      "NewQuery",      "Explain",
                                                             "Profile",  "CandidatePhase", "ExecutionPhase"};
  static inline const std::string& getTypeName(Type type) { return type_names[type]; }

  explicit ServerEvent(Type event_type) : type(event_type) {}

  friend inline std::ostream& operator<<(std::ostream& os, ServerEvent& event) {
    os << getTypeName(event.type);
    return os;
  }
};

enum class CandidatePruningStrategy : uint16_t { None = 0, Adaptive, LDF, NLF, CFL, DAF, GQL, TSO };
enum class CompressionStrategy : uint16_t { None = 0, Static, Dynamic };
enum class QueryMode : uint16_t { Execute = 0, Profile, Explain };

// TODO(tatiana): this is workaround as we do not implement the order now
inline std::vector<QueryVertexID> getOrder(const std::string& order_str, uint32_t size) {
  std::vector<QueryVertexID> order;
  if (order_str.empty()) return order;
  order.reserve(size);
  QueryVertexID v = 0;
  for (auto c : order_str) {
    if (c == ' ') {
      order.push_back(v);
      v = 0;
    } else {
      v = v * 10 + (c - '0');
    }
  }
  order.push_back(v);
  CHECK_EQ(order.size(), size);
  return order;
}

class QueryConfig {
 public:
  static inline const unordered_map<std::string, CompressionStrategy> compression_strategies = {
      {"none", CompressionStrategy::None},
      {"static", CompressionStrategy::Static},
      {"dynamic", CompressionStrategy::Dynamic}};
  static inline const unordered_map<std::string, CandidatePruningStrategy> candidate_pruning_strategies = {
      {"none", CandidatePruningStrategy::None}, {"adaptive", CandidatePruningStrategy::Adaptive},
      {"ldf", CandidatePruningStrategy::LDF},   {"nlf", CandidatePruningStrategy::NLF},
      {"cfl", CandidatePruningStrategy::CFL},   {"daf", CandidatePruningStrategy::DAF},
      {"gql", CandidatePruningStrategy::GQL}};
  static inline const unordered_map<std::string, QueryMode> modes = {
      {"execute", QueryMode::Execute}, {"profile", QueryMode::Profile}, {"explain", QueryMode::Explain}};

  std::string matching_order;
  CandidatePruningStrategy candidate_pruning_strategy = CandidatePruningStrategy::CFL;
  CompressionStrategy compression_strategy = CompressionStrategy::Dynamic;
  bool use_auxiliary_index = false;
  bool use_partitioned_graph = true;
  std::string output = "count";
  uint64_t limit = ~0ull;
  QueryMode mode = QueryMode::Execute;

  explicit QueryConfig(const std::string& config_str = "") {
    uint32_t tok_start = 0;
    for (uint32_t i = 0; i < config_str.size();) {
      while (config_str[i] != '=' && i < config_str.size()) {
        ++i;
      }
      if (i == config_str.size()) {
        throw std::runtime_error("Query config parse error. Expected = in " + config_str.substr(tok_start) +
                                 ". Query config format: key=value[,key=value]...");
      }
      std::string_view key(config_str.data() + tok_start, i - tok_start);
      tok_start = i + 1;
      while (config_str[i] != ',' && i < config_str.size()) {
        ++i;
      }
      std::string value(config_str.data() + tok_start, i - tok_start);
      tok_start = i + 1;
      if (key == "cps" || key == "candidate_pruning_strategy") {
        validateConfig(candidate_pruning_strategy, value, candidate_pruning_strategies, "candidate pruning strategy");
      } else if (key == "mo" || key == "matching_order") {
        matching_order = value;
        // TODO(tatiana): use strategy name instead of actual matching order
      } else if (key == "cs" || key == "compression_strategy") {
        validateConfig(compression_strategy, value, compression_strategies, "compression strategy");
      } else if (key == "mode") {
        validateConfig(mode, value, modes, "query mode");
      } else if (key == "limit") {
        limit = std::stoull(value);
      } else if (key == "use_partitioned_graph" || key == "upg") {
        use_partitioned_graph = value == "true" || value == "1";
      } else if (key == "use_auxiliary_index" || key == "uai") {
        use_auxiliary_index = value == "true" || value == "1";
      } else if (key == "output") {
        output = value;
      }
    }
  }

 private:
  template <typename T, typename Container>
  void validateConfig(T& conf, const std::string& val_str, const Container& options, const std::string& name) {
    auto pos = options.find(val_str);
    if (pos == options.end()) {
      std::stringstream ss;
      ss << "Invalid " << name << ' ' << val_str << ". {";
      for (auto& val : options) {
        ss << val.first << ",";
      }
      ss << "}";
      throw std::runtime_error(ss.str());
    } else {
      conf = pos->second;
    }
  }
};

struct QueryContext {
  QueryGraph query_graph;
  QueryConfig query_config;
  GraphBase* data_graph;
  GraphMetadata* graph_metadata;

  QueryContext(QueryGraph&& q, QueryConfig&& config, GraphBase* g, GraphMetadata* metadata)
      : query_graph(std::move(q)), query_config(std::move(config)), data_graph(g), graph_metadata(metadata) {}

  void operator=(QueryContext&& context) {
    query_graph = std::move(context.query_graph);
    query_config = std::move(context.query_config);
    data_graph = context.data_graph;
    graph_metadata = context.graph_metadata;
  }
};

struct QueryResult {
  double filter_time = 0;
  double plan_time = 0;
  double enumerate_time = 0;
  double elapsed_execution_time = 0;
  uint64_t embedding_count = 0;
  std::string matching_order;

  inline std::string toString() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }

  friend std::ostream& operator<<(std::ostream& out, const QueryResult& res) {
    return out << res.elapsed_execution_time << ',' << res.filter_time << ',' << res.plan_time << ','
               << res.enumerate_time << ',' << res.embedding_count << ',' << res.matching_order;
  }
};

struct ProfileInfo {
  uint64_t total_input_size = 0;
  uint64_t total_output_size = 0;
  uint64_t total_num_input_subgraphs = 0;
  uint64_t total_num_output_subgraphs = 0;
  double total_time_in_milliseconds = 0;
  uint64_t intersection_count = 0;
  uint64_t total_intersection_input_size = 0;
  uint64_t total_intersection_output_size = 0;
  /** The minimal number of intersection needed if all intersection function call results can be cached */
  uint64_t distinct_intersection_count = 0;

  void updateIntersectInfo(uint32_t input_size, uint32_t output_size) {
    ++intersection_count;
    total_intersection_input_size += input_size;
    total_intersection_output_size += output_size;
  }
};

}  // namespace circinus

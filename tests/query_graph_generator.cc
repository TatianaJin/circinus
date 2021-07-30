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
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <vector>
#include "gflags/gflags.h"
#include "graph/graph.h"
#include "graph/types.h"

DEFINE_uint64(qg_cnt, 200, "target query graph count");
DEFINE_uint64(v_cnt, 4, "vertex count in each query");
DEFINE_bool(if_dense, false, "if dense query graph");
DEFINE_string(output_dir, "/data/share/users/qlma/query-graph-output/", "where to put the result");
DEFINE_string(data_graph, "/data/share/project/haxe/data/subgraph_matching_datasets/yeast/data_graph/yeast.graph",
              "input data graph");
namespace circinus {

class QueryGraphGenerator {
 private:
  Graph g_;
  uint32_t target_query_graph_cnt_;
  uint32_t target_vertex_cnt_;
  bool if_dense_;
  std::string output_dir_;
  uint32_t max_random_walk_step;
  uint32_t attempt_cnt_;
  uint32_t success_cnt_;
  uint32_t attempt_cnt_bound_;
  std::set<VertexID> vset_;
  std::set<std::pair<VertexID, VertexID>> eset_;
  double avg_deg_;
  
  const uint32_t dense_sparse_option_bound_ = 5;
  // there is an option for dense/sparse only when the amount of query vertices reach this bound
  const uint32_t dense_deg_bound_ = 3;
  // the minimum average degree of a dense query graph

 public:
  QueryGraphGenerator(std::string data_graph, uint32_t target_query_graph_cnt, uint32_t target_vertex_cnt,
                      bool if_dense, std::string output_dir)
      : g_(data_graph),
        target_query_graph_cnt_(target_query_graph_cnt),
        target_vertex_cnt_(target_vertex_cnt),
        if_dense_(if_dense),
        output_dir_(output_dir),
        max_random_walk_step(5 * target_vertex_cnt),
        attempt_cnt_(0),
        success_cnt_(0),
        attempt_cnt_bound_(1e8),
        avg_deg_(0) {}
  void generate() {
    std::mt19937 mt_rand(std::chrono::system_clock::now().time_since_epoch().count());
    vset_.clear();
    eset_.clear();
    uint32_t step_cnt = 0;
    VertexID current_vertex = mt_rand() % g_.getNumVertices();
    vset_.insert(current_vertex);
    std::set<std::pair<VertexID, VertexID>> extra_edge_set;
    while (step_cnt < max_random_walk_step) {
      auto[neighbors, neighbors_cnt] = g_.getOutNeighbors(current_vertex);
      if (neighbors_cnt == 0) break;
      for (auto neighbor_v : vset_) {
        if (std::binary_search(neighbors, neighbors + neighbors_cnt, neighbor_v)) {
          if (current_vertex < neighbor_v)
            extra_edge_set.insert({current_vertex, neighbor_v});
          else
            extra_edge_set.insert({neighbor_v, current_vertex});
        }
      }
      if (vset_.size() >= target_vertex_cnt_) break;
      VertexID new_vertex = neighbors[mt_rand() % neighbors_cnt];
      vset_.insert(new_vertex);
      if (current_vertex < new_vertex)
        eset_.insert({current_vertex, new_vertex});
      else
        eset_.insert({new_vertex, current_vertex});
      current_vertex = new_vertex;
      ++step_cnt;
    }
    std::vector<std::pair<VertexID, VertexID>> diff_vec;
    std::set_difference(extra_edge_set.begin(), extra_edge_set.end(), eset_.begin(), eset_.end(),
                        std::back_inserter(diff_vec));
    int all_edge_cnt = diff_vec.size() + eset_.size();
    avg_deg_ += 2.0 * all_edge_cnt / vset_.size();
    // calculate the maximum possible degree (if all edges between pairs in vset_ were selected)
    if (vset_.size() < target_vertex_cnt_ || !if_dense_ || target_vertex_cnt_ < dense_sparse_option_bound_) return;
    // can be changed according to generating rules
    // deal with the situation: 1.enough vertex 2.target dense qg(normally means not enough edges)
    int d = 2 * all_edge_cnt - dense_deg_bound_ * target_vertex_cnt_;
    if (d < 0) return;
    d /= 2;
    std::random_shuffle(diff_vec.begin(), diff_vec.end());
    int new_edge_cnt_bound = diff_vec.size();
    new_edge_cnt_bound -= mt_rand() % (d + 1);
    for (int i = 0; i < new_edge_cnt_bound; ++i) {
      eset_.insert(diff_vec[i]);
    }
  }
  bool check() {
    // can be changed according to generating rules
    if (vset_.size() != target_vertex_cnt_) return 0;
    if (target_vertex_cnt_ < dense_sparse_option_bound_) return 1;
    if (if_dense_) {
      if (2 * eset_.size() < dense_deg_bound_ * target_vertex_cnt_) return 0;
    } else {
      if (2 * eset_.size() >= dense_deg_bound_ * target_vertex_cnt_) return 0;
    }
    return 1;
  }
  std::string getFilename() {
    // can be changed according to generating rules
    std::stringstream output_filename;
    output_filename << output_dir_ << "/query_";
    if (target_vertex_cnt_ < dense_sparse_option_bound_ || if_dense_)
      output_filename << "dense_";
    else
      output_filename << "sparse_";
    output_filename << std::to_string(target_vertex_cnt_) << "_";
    output_filename << std::to_string(success_cnt_) << ".graph";
    return output_filename.str();
  }
  void save() {
    std::string fn = getFilename();
    std::map<VertexID, int> id_map;
    std::map<int, LabelID> label_map;
    std::map<int, int> degree_map;
    std::vector<std::pair<int, int>> elist;
    int newid = 0;
    for (auto& v : vset_) {
      id_map[v] = newid;
      label_map[newid] = g_.getVertexLabel(v);
      ++newid;
    }
    for (auto & [ v1, v2 ] : eset_) {
      int id1 = id_map[v1], id2 = id_map[v2];
      ++degree_map[id1];
      ++degree_map[id2];
      elist.push_back({id1, id2});
    }
    std::ofstream out(fn);
    out << "t " << target_vertex_cnt_ << " " << elist.size() << "\n";
    for (int i = 0; i < newid; ++i) {
      out << "v " << i << " " << label_map[i] << " " << degree_map[i] << "\n";
    }
    for (auto & [ v1, v2 ] : elist) {
      out << "e " << v1 << " " << v2 << " 0\n";
    }
    out.close();
  }
  void run() {
    while (success_cnt_ < target_query_graph_cnt_ && attempt_cnt_ < attempt_cnt_bound_) {
      generate();
      ++attempt_cnt_;
      if (check()) {
        ++success_cnt_;
        save();
      }
    }
    std::cout << "Tried " << attempt_cnt_ << " times to get " << success_cnt_ << " queries.\n";
    std::cout << "Avg degree of all tried random walk induced subgraphs: " << avg_deg_ / attempt_cnt_ << "\n";
  }
};
}  // namespace circinus

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, false);
  circinus::QueryGraphGenerator qgg(FLAGS_data_graph, FLAGS_qg_cnt, FLAGS_v_cnt, FLAGS_if_dense, FLAGS_output_dir);
  qgg.run();
}

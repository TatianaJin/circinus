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
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include "graph/graph.h"
#include "graph/types.h"
#include "gflags/gflags.h"

DEFINE_uint64(qg_cnt, 200, "target query graph count");
DEFINE_uint64(v_cnt, 4, "vertex count in each query");
DEFINE_bool(if_dense, false, "if dense query graph");

namespace circinus {

class QueryGraphGenerator{
 private:
  Graph g_;
  uint32_t target_query_graph_cnt_;
  uint32_t target_vertex_cnt_;
  bool if_dense_;
  uint32_t max_random_walk_step;
  uint32_t trial_cnt_;
  uint32_t success_cnt_;
	uint32_t trial_cnt_bound_;
  std::set<VertexID> vset_;
  std::set< std::pair<VertexID,VertexID> > eset_;
 public:
  QueryGraphGenerator(Graph g,uint32_t target_query_graph_cnt,uint32_t target_vertex_cnt,bool if_dense):
    g_(g),
    target_query_graph_cnt_(target_query_graph_cnt),
    target_vertex_cnt_(target_vertex_cnt),
    if_dense_(if_dense),
    max_random_walk_step(5*target_vertex_cnt),
    trial_cnt_(0),
    success_cnt_(0),
		trial_cnt_bound_(1e8)
    {}
  void generate()
  {
    std::mt19937 mt_rand(std::random_device{}());
    vset_.clear();
    eset_.clear();
    uint32_t step_cnt=0;
    VertexID current_vertex=mt_rand()%g_.getNumVertices();
		vset_.insert(current_vertex);
    while(step_cnt<max_random_walk_step&&vset_.size()<target_vertex_cnt_)
    {
      auto [neighbors,neighbors_cnt]=g_.getOutNeighbors(current_vertex);
			if(neighbors_cnt==0)break;
      VertexID new_vertex = neighbors[mt_rand()%neighbors_cnt];
      vset_.insert(new_vertex);
      eset_.insert({current_vertex,new_vertex});
      eset_.insert({new_vertex,current_vertex});
      current_vertex=new_vertex;
      ++step_cnt;
    }
    if(vset_.size()<target_vertex_cnt_||!if_dense_||target_vertex_cnt_<8)return;
    //deal with the situation: 1.enough vertex 2.target dense qg(normally means not enough edges)
    std::vector<VertexID> permutation1(vset_.begin(),vset_.end()),permutation2(vset_.begin(),vset_.end());
    std::random_shuffle(permutation1.begin(),permutation1.end());
    std::random_shuffle(permutation2.begin(),permutation2.end());
    uint32_t try_new_edge_cnt=0;
    uint32_t try_new_edge_cnt_bound=0; 
    for(VertexID start_v:permutation1)
    {
      for(VertexID end_v:permutation2)
      {
        ++try_new_edge_cnt;
        if(try_new_edge_cnt_bound&&try_new_edge_cnt>try_new_edge_cnt_bound)return;
        auto [neighbors,neighbors_cnt]=g_.getOutNeighbors(start_v);
        if(!neighbors_cnt)continue;
        int l=0,r=neighbors_cnt-1;
        while(l<r)
        {
          int mid=(l+r)/2;
          if(neighbors[mid]>=end_v)r=mid;
          else l=mid+1;
        }
        if(neighbors[l]==end_v)
        {
          eset_.find({start_v,end_v});
          eset_.find({end_v,start_v});
          if(eset_.size()/vset_.size()>=3&&!try_new_edge_cnt_bound)
          {
            try_new_edge_cnt_bound=try_new_edge_cnt+mt_rand()%(target_vertex_cnt_*target_vertex_cnt_-try_new_edge_cnt+1);
          }
        }
      }
    }
  }
  bool check()
  {
      if(vset_.size()!=target_vertex_cnt_)return 0;
			if(target_vertex_cnt_<8)return 1;
      if(if_dense_)
      {
        if(eset_.size()/vset_.size()<3)return 0;
      }
      else
      {
        if(eset_.size()/vset_.size()>=3)return 0;
      }
			
      return 1;
  }
  void save(){}
  void run()
  {
    while(success_cnt_<target_query_graph_cnt_&&trial_cnt_<trial_cnt_bound_)
    {
      generate();
      ++trial_cnt_;
      if(check())
      {
        ++success_cnt_;
        save();
      }
    }
    std::cout<<"Tried "<<trial_cnt_<<" times to get "<<success_cnt_<<" queries.\n";
  }
};

}

const char data_graph_dir[] = "/data/share/project/haxe/data/subgraph_matching_datasets/yeast/data_graph/yeast.graph";
int main(int argc, char** argv)
{
  gflags::ParseCommandLineFlags(&argc, &argv, false);
  auto data_graph_dir_str = std::string(data_graph_dir);
  circinus::Graph g(data_graph_dir_str);
  circinus::QueryGraphGenerator qgg(g,FLAGS_qg_cnt,FLAGS_v_cnt,FLAGS_if_dense);
  qgg.run();
}

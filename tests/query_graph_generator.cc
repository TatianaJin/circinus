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
    std::set< std::pair<VertexID,VertexID> > extra_edge_set;
    while(step_cnt<max_random_walk_step)
    {
      auto [neighbors,neighbors_cnt]=g_.getOutNeighbors(current_vertex);
			if(neighbors_cnt==0)break;
      if(if_dense_)
      {
        for(auto neighbor_v : vset_)
        {
          int l=0,r=neighbors_cnt-1;
          while(l<r)
          {
            int mid=(l+r)/2;
            if(neighbors[mid]>=end_v)r=mid;
            else l=mid+1;
          }
          if(neighbors[l]==neighbor_v)
          {
            extra_edge_set.insert({current_vertex,neighbor_v});
            extra_edge_set.insert({neighbor_v,current_vertex});
          }
        }
      }
      if(vset_.size()>=target_vertex_cnt_)break;
      VertexID new_vertex = neighbors[mt_rand()%neighbors_cnt];
      vset_.insert(new_vertex);
      eset_.insert({current_vertex,new_vertex});
      eset_.insert({new_vertex,current_vertex});
      current_vertex=new_vertex;
      ++step_cnt;
    }
    if(vset_.size()<target_vertex_cnt_||!if_dense_||target_vertex_cnt_<8)return;
    //deal with the situation: 1.enough vertex 2.target dense qg(normally means not enough edges)
    std::vector< std::pair<VertexID,VertexID> > diff_vec;
    std::set_difference(extra_edge_set.begin(),extra_edge_set.end(),eset_.begin(),eset_.end(),std::back_inserter(diff_vec));
    uint32_t all_edge_cnt=diff_vec.size()+eset_.size();
    uint32_t d=all_edge_cnt-3*target_vertex_cnt_;
    if(d<0)return;
    d/=2;
    std::random_shuffle(diff_vec.begin(),diff_vec.end());
    uint32_t new_edge_cnt=diff_vec.size()/2;
    if(d>0)new_edge_cnt-=mt_rand()%(d+1);
    for(uint32_t i=0;i<new_edge_cnt;++i)
    {
      auto p = diff_vec[i];
      if(eset_.find(p)==eset_.end())
      {
        eset_.insert(p);
        eset_.insert({p.second,p.first});
      }
      else --i;
    }
  }
  bool check()
  {
      if(vset_.size()!=target_vertex_cnt_)return 0;
			if(target_vertex_cnt_<8)return 1;
      if(if_dense_)
      {
        if(eset_.size()<3*target_vertex_cnt_)return 0;
      }
      else
      {
        if(eset_.size()>=3*target_vertex_cnt_)return 0;
      }
      return 1;
  }
  void save(){}
  void run()
  {
    double avg_deg=0;
    while(success_cnt_<target_query_graph_cnt_&&trial_cnt_<trial_cnt_bound_)
    {
      generate();
      if(vset_.size())avg_deg+=1.0*eset_.size()/vset_.size();
      ++trial_cnt_;
      if(check())
      {
        ++success_cnt_;
        save();
      }
    }
    std::cout<<"Tried "<<trial_cnt_<<" times to get "<<success_cnt_<<" queries.\n";
    std::cout<<"Avg degree of all tried subgraph: "<<avg_deg<<"\n";
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

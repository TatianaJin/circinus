#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"
#ifdef WITH_GPERF
#include "gperftools/profiler.h"
#endif
#include "gtest/gtest.h"

#include "exec/thread_pool.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/query_graph.h"
#include "graph/bipartite_graph.h"
#include "ops/filters.h"
#include "ops/logical_filters.h"
#include "ops/operators.h"
#include "ops/order.h"
#include "ops/scans.h"
#include "ops/types.h"
#include "plan/execution_plan.h"
#include "plan/naive_planner.h"
#include "utils/flags.h"
#include "utils/hashmap.h"
#include "utils/profiler.h"
#include "utils/utils.h"

using circinus::CompressedSubgraphs;
using circinus::ExecutionConfig;
using circinus::ExecutionPlan;
using circinus::Graph;
using circinus::BipartiteGraph;
using circinus::GraphType;
using circinus::GraphMetadata;
using circinus::NaivePlanner;
using circinus::QueryGraph;
using circinus::QueryVertexID;
using circinus::Task;
using circinus::ThreadPool;
using circinus::VertexID;
using circinus::Profiler;
using circinus::CoverNode;
using circinus::QueryType;
using circinus::TraverseOperator;
using circinus::INVALID_VERTEX_ID;

// logical filter
using circinus::LogicalCFLFilter;
using circinus::LogicalGQLFilter;
using circinus::LogicalNLFFilter;
using circinus::LogicalTSOFilter;
using circinus::LogicalDPISOFilter;
using circinus::LogicalNeighborhoodFilter;

// physical filter
using circinus::NeighborhoodFilter;
using circinus::NLFFilter;
using circinus::GQLFilter;
namespace circinus{
class FilterAndOrder
{
  private:
    std::unique_ptr<LogicalNeighborhoodFilter> logical_filter_;
    std::map<std::pair<QueryVertexID,QueryVertexID>,BipartiteGraph> bg_map_;
    bool filtered=0;
    std::vector<std::vector<VertexID>> candidates_sets_;
    const Graph* g_pointer_;
    const QueryGraph* q_pointer_;
    std::string filter_string_;
  public:
    FilterAndOrder(const Graph* g_pointer,const QueryGraph* q_pointer,std::string filter_string):
                   g_pointer_(g_pointer),q_pointer_(q_pointer),filter_string_(filter_string){}

    const BipartiteGraph* getBipartiteGraph(QueryVertexID v1,QueryVertexID v2)
    {
      assert(filtered);
      std::pair<QueryVertexID,QueryVertexID> p(v1,v2);
      auto bg = bg_map_.find(p);
      if(bg==bg_map_.end())
      {
        BipartiteGraph newbg(v1,v2);
        newbg.populateGraph(g_pointer_,&candidates_sets_);
        auto res=bg_map_.insert({p,std::move(newbg)});
        bg=res.first;
      }
      return &(bg->second);
    }

    std::vector<std::vector<VertexID>> getCandidateSets() { 
      if(filtered)return candidates_sets_;
      const Graph g=*g_pointer_;
      const QueryGraph q=*q_pointer_;
      candidates_sets_.resize(q.getNumVertices());
      std::vector<VertexID> candidate_size(q.getNumVertices());
      ExecutionConfig config;
      for (uint32_t v = 0; v < q.getNumVertices(); ++v) {
        config.setInputSize(g.getVertexCardinalityByLabel(q.getVertexLabel(v)));
        auto scan = circinus::Scan::newLDFScan(q.getVertexLabel(v), q.getVertexOutDegree(v), 0, config, 1);
        if (filter_string_.compare("ldf") && filter_string_.compare("dpiso")) {
          scan->addFilter(std::make_unique<NLFFilter>(&q, v));
        }
        auto scan_ctx = scan->initScanContext(0);
        scan->scan(&g, &scan_ctx);
        candidates_sets_[v] = std::move(scan_ctx.candidates);
        candidate_size[v] = candidates_sets_[v].size();
      }

      if (filter_string_.compare("ldf") && filter_string_.compare("nlf")) {
        auto metadata = GraphMetadata(g);
        if (filter_string_.compare("cfl")==0) {
          logical_filter_ = std::make_unique<LogicalCFLFilter>(metadata, &q, candidate_size);
        } else if (filter_string_.compare("dpiso")==0) {
          logical_filter_ = std::make_unique<LogicalDPISOFilter>(metadata, &q, candidate_size);
        } else if (filter_string_.compare("tso")==0) {
          logical_filter_ = std::make_unique<LogicalTSOFilter>(metadata, &q, candidate_size);
        } else if (filter_string_.compare("gql")==0) {
          logical_filter_ = std::make_unique<LogicalGQLFilter>(&q);
        }
        auto physical_filters = logical_filter_->toPhysicalOperators(metadata, config);
        for (auto& filter : physical_filters) {
          QueryVertexID query_vertex = filter->getQueryVertex();
          filter->setInputSize(candidate_size[query_vertex]);
          auto filter_ctx = filter->initFilterContext(0);
          filter->filter(&g, &candidates_sets_, &filter_ctx);
          candidates_sets_[query_vertex].erase(std::remove_if(candidates_sets_[query_vertex].begin(), candidates_sets_[query_vertex].end(),
                                                        [invalid_vertex_id = INVALID_VERTEX_ID](VertexID & candidate) {
                                                          return candidate == invalid_vertex_id;
                                                        }),
                                        candidates_sets_[query_vertex].end());
          candidate_size[query_vertex] = candidates_sets_[query_vertex].size();
        }
      }
      filtered=1;
      return candidates_sets_;
    }
  QueryVertexID selectGQLStartVertex(const QueryGraph *query_graph)
  {
    QueryVertexID start_vertex=0;
    QueryVertexID qg_v_cnt= query_graph->getNumVertices();
    for(QueryVertexID i=1;i<qg_v_cnt;++i)
    {
      QueryVertexID cur_vertex=i;
      size_t size1=candidates_sets_[cur_vertex].size(),size2=candidates_sets_[start_vertex].size();
      if(size1<size2)
      {
        start_vertex = cur_vertex;
      }
      else
      if (size1==size2&&query_graph->getVertexOutDegree(cur_vertex)>query_graph->getVertexOutDegree(start_vertex))
      {
        start_vertex = cur_vertex;
      }
    }
    return start_vertex;
  }

  void updateValidVertices(QueryVertexID query_vertex,std::vector<bool> &visited,std::vector<bool> &adjacent)
  {
    visited[query_vertex] = true;
    auto [nbrs,cnt]=q_pointer_->getOutNeighbors(query_vertex);
    for(uint32_t i=0;i<cnt;++i)
    {
      adjacent[nbrs[i]]=true;
    }
  }

  std::vector<QueryVertexID> getGQLOrder(
    const Graph *data_graph, const QueryGraph *query_graph
  )
  {
    QueryVertexID qg_v_cnt= query_graph->getNumVertices();
    std::vector<bool> visited_vertices(qg_v_cnt, false);
    std::vector<bool> adjacent_vertices(qg_v_cnt, false);
    std::vector<QueryVertexID> order(qg_v_cnt);

    QueryVertexID start_vertex= selectGQLStartVertex(query_graph);
    order[0]=start_vertex;
    updateValidVertices(start_vertex,visited_vertices,adjacent_vertices);

    for(QueryVertexID i=1;i<qg_v_cnt;++i)
    {
      QueryVertexID next_vertex;
      QueryVertexID min_value = data_graph->getNumVertices()+1;
      for(QueryVertexID j=0;j<qg_v_cnt;++j)
      {
        QueryVertexID cur_vertex=j;
        if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
          size_t cnt=candidates_sets_[cur_vertex].size();
          if(cnt<min_value)
          {
            min_value=cnt;
            next_vertex = cur_vertex;
          }
          else if (cnt==min_value&&query_graph->getVertexOutDegree(cur_vertex)>query_graph->getVertexOutDegree(next_vertex))
          {
            next_vertex = cur_vertex;
          }
        }
      }
      updateValidVertices(next_vertex,visited_vertices,adjacent_vertices);
      order[i]=next_vertex;
    }
    return order;
  }
  
  std::vector<QueryVertexID> getDPISOOrder(
    const QueryGraph *query_graph,
    std::vector<QueryVertexID> bfs_order
  )
  {
    QueryVertexID qg_v_cnt= query_graph->getNumVertices();
    std::vector<QueryVertexID> order(qg_v_cnt);
    for(QueryVertexID i=0;i<qg_v_cnt;++i)
    {
      order[i]=bfs_order[i];
    }
    return order;
  }

  void estimatePathEmbeddsingsNum(
    std::vector<QueryVertexID> &path,
    std::vector<size_t> &estimated_embeddings_num)
    {
      assert(path.size() > 1);
      std::vector<size_t> parent;
      std::vector<size_t> children;

      size_t begin=path.size() - 2,end=path.size() - 1;

      estimated_embeddings_num.resize(path.size() - 1);
      auto last_edge =getBipartiteGraph(path[begin],path[end]);
      children.resize(last_edge->getNumVertices());
      
      size_t sum = 0;
      for (auto& v : candidates_sets_[path[begin]]) {
          int offset=last_edge->getOffset(v);
          children[offset] = last_edge->getVertexOutDegree(v);
          sum += children[offset];
      }

      estimated_embeddings_num[begin] = sum;

      for (int i = begin; i >= 1; --i) {
        begin = path[i - 1];
        end = path[i];
        auto edge = getBipartiteGraph(path[begin],path[end]);
        parent.resize(edge->getNumVertices());

        sum=0;
        for(auto& v : candidates_sets_[path[begin]])
        {
          size_t local_sum = 0;
          auto& [nbrs,cnt]=edge.getOutNeighbors(v);
          for(uint32_t j=0;j<cnt;++j)
          {
            auto nbr=nbrs[j];
            local_sum += children[last_edge->getOffset(nbr)];
          }
          parent[edge->getOffset(v)]=local_sum;
          sum+=local_sum;
        }

        estimated_embeddings_num[i - 1] = sum;
        parent.swap(children);

        last_edge=edge;
      }
    }
  QueryVertexID generateNoneTreeEdgesCount(
    const QueryGraph *query_graph,
    const std::vector<TreeNode>& tree_node,
    std::vector<QueryVertexID> &path)
    {
      auto non_tree_edge_count = query_graph->getVertexOutDegree(path[0]) - tree_node[path[0]].children_.size();
      for (size_t i = 1; i < path.size(); ++i) {
        auto vertex = path[i];
        non_tree_edge_count += query_graph->getVertexOutDegree(vertex) - tree_node[vertex].children_.size() - 1;
      }

      return non_tree_edge_count;
    }

  void generateRootToLeafPaths(
    const std::vector<TreeNode>& tree_node,
    QueryVertexID cur_vertex,
    std::vector<QueryVertexID> &cur_path,
    std::vector<std::vector<QueryVertexID>>&paths )
  {
    TreeNode& cur_node = tree_node[cur_vertex];
    cur_path.push_back(cur_vertex);
    if(cur_node.children_.size()==0)
    {
      paths.emplace_back(cur_path);
    }
    else
    {
      for(auto child:cur_node.children_)
      {
        generateRootToLeafPaths(tree_node, next_vertex, cur_path, paths);
      }
    }
    cur_path.pop_back();
  }
  std::vector<QueryVertexID> getTSOOrder(
    const QueryGraph *query_graph,
    std::map<std::pair<QueryVertexID,QueryVertexID>,BipartiteGraph> bg_map,
    const std::vector<TreeNode>& tree,
    std::vector<QueryVertexID> dfs_order
    )
  {
    QueryVertexID qg_v_cnt= query_graph->getNumVertices();
    std::vector<std::vector<QueryVertexID>>paths;
    paths.reserve(qg_v_cnt);

    std::vector<QueryVertexID> single_path;
    single_path.reserve(qg_v_cnt);

    generateRootToLeafPaths(tree,dfs_order[0],single_path,paths);
    std::vector<std::pair<double, std::vector<QueryVertexID>*>> path_orders;
    for(auto& path:paths)
    {
      std::vector<size_t> estimated_embeddings_num;
      VertexID non_tree_edges_count = generateNoneTreeEdgesCount(query_graph,tree,path);
      estimatePathEmbeddsingsNum(path, estimated_embeddings_num);
      double score = estimated_embeddings_num[0] / (double) (non_tree_edges_count + 1);
      path_orders.emplace_back(std::make_pair(score, &path));
    }
    std::sort(path_orders.begin(), path_orders.end());
    std::vector<bool> visited_vertices(qg_v_cnt, false);
    std::vector<QueryVertexID> order;
    order.reserve(qg_v_cnt);
    for (auto& path : path_orders) {
      for(auto v:*(path.second))
      {
        if(!visited_vertices[v])
        {
          order.push_back(v);
          visited_vertices[v]=true;
        }
      }
    return order;
  }
  
};

}


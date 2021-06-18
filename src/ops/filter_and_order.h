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
    std::string filter_string_;
  public:
    FilterAndOrder(std::string filter_string):filter_string_(filter_string){}
    std::vector<std::vector<VertexID>> getCandidateSets(const Graph& g, const QueryGraph& q) { 
      std::vector<std::vector<VertexID>> candidates(q.getNumVertices());
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
        candidates[v] = std::move(scan_ctx.candidates);
        candidate_size[v] = candidates[v].size();
        // LOG(INFO) << "query vertex " << v << ' ' << scan->toString();
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
          filter->filter(&g, &candidates, &filter_ctx);
          candidates[query_vertex].erase(std::remove_if(candidates[query_vertex].begin(), candidates[query_vertex].end(),
                                                        [invalid_vertex_id = INVALID_VERTEX_ID](VertexID & candidate) {
                                                          return candidate == invalid_vertex_id;
                                                        }),
                                        candidates[query_vertex].end());
          candidate_size[query_vertex] = candidates[query_vertex].size();
        }
      }
      return candidates;
    }
};

}


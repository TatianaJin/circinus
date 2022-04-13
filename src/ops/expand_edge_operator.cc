#include "ops/expand_edge_operator.h"

#include <algorithm>
#include <chrono>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/intersect.h"
#include "graph/compressed_subgraphs.h"
#include "graph/graph.h"
#include "graph/graph_view.h"
#include "graph/types.h"
#include "graph/vertex_set.h"
#include "ops/filters/subgraph_filter.h"
#include "ops/traverse_operator.h"
#include "ops/traverse_operator_utils.h"
#include "ops/types.h"
#include "utils/hashmap.h"

namespace circinus {

#define CONSTRUCT(name)                                                                                                \
  ExpandEdge##name##Operator(uint32_t parent_index, uint32_t target_index, QueryVertexID parent, QueryVertexID target, \
                             const std::vector<uint32_t>& same_label_key_indices,                                      \
                             const std::vector<uint32_t>& same_label_set_indices, uint64_t set_pruning_threshold,      \
                             std::unique_ptr<SubgraphFilter>&& filter)                                                 \
      : ExpandEdgeOperator(parent_index, target_index, parent, target, same_label_key_indices, same_label_set_indices, \
                           set_pruning_threshold, std::move(filter)) {}

template <typename G, bool intersect_candidates>
class ExpandEdgeKeyToSetOperator : public ExpandEdgeOperator {
 public:
  CONSTRUCT(KeyToSet)

  uint32_t expand(uint32_t cap, TraverseContext* ctx) const override {
    uint32_t n = 0;
    for (; n < cap && ctx->hasNextInput(); ctx->nextInput()) {
      n += expandInner<QueryType::Execute>(ctx->getCurrentInput(), ctx);
    }
    return n;
  }

  uint32_t expandAndProfileInner(uint32_t cap, TraverseContext* ctx) const override {
    uint32_t old_input_index = ctx->getInputIndex();
    uint32_t n = 0;
    for (; n < cap && ctx->hasNextInput(); ctx->nextInput()) {
      n += (ctx->getQueryType() == QueryType::Profile)
               ? expandInner<QueryType::Profile>(ctx->getCurrentInput(), ctx)
               : (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect
                      ? expandInner<QueryType::ProfileCandidateSIEffect>(ctx->getCurrentInput(), ctx)
                      : expandInner<QueryType::ProfileWithMiniIntersection>(ctx->getCurrentInput(), ctx));
      ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
    }
    ctx->intersection_count += ctx->getInputIndex() - old_input_index;
    return n;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandEdgeKeyToSetOperator";
    toStringInner(ss);
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first, input_key_size.second + 1};
  }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    return std::make_unique<ExpandEdgeTraverseContext>(candidates, graph, outputs, profile, intersect_candidates,
                                                       candidate_hashmap);
  }

 private:
  /** @returns True if one CompressedSubgraphs is generated, else false. */
  template <QueryType profile>
  inline bool expandInner(const CompressedSubgraphs& input, TraverseContext* ctx) const {
    VertexSet target_set;
    if (canReuseSet()) {
      target_set = input.getSet(reusable_set_index_);
      if (target_set->size() == 1) return false;
      filterTargets(*target_set, input);
      if (target_set->empty()) return false;
    } else {
      std::vector<VertexID> targets;
      auto parent_match = input.getKeyVal(parent_index_);
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      expandFromParent<G, profile, intersect_candidates>(ctx, parent_match, exceptions, &targets, input);

      if (targets.empty()) {
        return false;
      }
      target_set = newVertexSet(targets);
    }
#ifdef USE_FILTER
    auto output = ctx->newOutput(input, std::move(target_set));
    if (filter(*output)) {  // actively prune existing sets
      ctx->popOutput();
      return false;
    }
#else
    auto output = ctx->newOutput(input, std::move(target_set), same_label_set_indices_, set_pruning_threshold_);
    if (output == nullptr) {  // actively prune existing sets
      return false;
    }
#endif
    return true;
  }
};  // class ExpandEdgeKeyToSetOperator

class ExpandEdgeKeyToKeyTraverseContext : public ExpandEdgeTraverseContext, public TargetBuffer {
 public:
  ExpandEdgeKeyToKeyTraverseContext(const CandidateSetView* candidates, const void* graph,
                                    std::vector<CompressedSubgraphs>* outputs, QueryType profile, bool hash_candidates,
                                    const unordered_set<VertexID>* candidate_hashmap)
      : ExpandEdgeTraverseContext(candidates, graph, outputs, profile, hash_candidates, candidate_hashmap) {}

  std::unique_ptr<TraverseContext> clone() const override {
    return std::make_unique<ExpandEdgeKeyToKeyTraverseContext>(*this);
  }
};

template <typename G, bool intersect_candidates>
class ExpandEdgeKeyToKeyOperator : public ExpandEdgeOperator {
 public:
  CONSTRUCT(KeyToKey)

  uint32_t expand(uint32_t cap, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(cap, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t cap, TraverseContext* ctx) const override {
    auto old_input_index = ctx->getInputIndex();
    uint32_t n = 0;
    if (ctx->getQueryType() == QueryType::Profile) {
      n = expandInner<QueryType::Profile>(cap, ctx);
    } else if (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect) {
      n = expandInner<QueryType::ProfileCandidateSIEffect>(cap, ctx);
    } else {
      CHECK(ctx->getQueryType() == QueryType::ProfileWithMiniIntersection) << "unknown query type "
                                                                           << (uint32_t)ctx->getQueryType();
      n = expandInner<QueryType::ProfileWithMiniIntersection>(cap, ctx);
    }
    ctx->intersection_count += (ctx->getInputIndex() - old_input_index) * (intersect_candidates);
    return n;
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandEdgeKeyToKeyOperator";
    toStringInner(ss);
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first + 1, input_key_size.second + 1};
  }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    return std::make_unique<ExpandEdgeKeyToKeyTraverseContext>(candidates, graph, outputs, profile,
                                                               intersect_candidates, candidate_hashmap);
  }

 private:
  template <QueryType profile>
  uint32_t expandInner(uint32_t cap, TraverseContext* base_ctx) const {
    uint32_t n = 0;
    auto ctx = dynamic_cast<ExpandEdgeKeyToKeyTraverseContext*>(base_ctx);
    while (true) {
      // if there are existing targets from the last input, consume first
      if (ctx->hasTarget()) {
        auto& input = ctx->getPreviousInput();
        while (ctx->hasTarget()) {
#ifdef USE_FILTER
          auto output =
              ctx->newOutput(input, ctx->currentTarget(), same_label_set_indices_, set_pruning_threshold_, false);
#else
          auto output = ctx->newOutput(input, ctx->currentTarget(), same_label_set_indices_, set_pruning_threshold_);
#endif
          ctx->nextTarget();
          if (output == nullptr) {
            continue;
          }
#ifdef USE_FILTER
          if (filter(*output)) {
            ctx->popOutput();
            continue;
          }
#endif
          if (++n == cap) {
            return n;
          }
        }
      }
      // return if all inputs in the current batch are consumed
      if (!ctx->hasNextInput()) {
        return n;
      }
      // consume the next input
      expandInner<profile>(ctx->getCurrentInput(), ctx);
      if (isProfileMode(profile)) {
        ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
      }
      ctx->nextInput();
    }
    return n;
  }

  template <QueryType profile>
  inline void expandInner(const CompressedSubgraphs& input, ExpandEdgeKeyToKeyTraverseContext* ctx) const {
    if (canReuseSet()) {
      auto target_view = *input.getSet(reusable_set_index_);
      if (target_view.size() == 1) {
        ctx->resetTargets();
      } else {
        filterTargets(target_view, input);
        ctx->setTargetView(target_view);
      }
    } else {
      auto& current_targets = ctx->resetTargets();
      auto parent_match = input.getKeyVal(parent_index_);
      auto exceptions = input.getExceptions(same_label_key_indices_, same_label_set_indices_);
      expandFromParent<G, profile, intersect_candidates>(ctx, parent_match, exceptions, &current_targets, input);
      ctx->resetTargetView();
    }
  }
};  // class ExpandEdgeKeyToKeyOperator

class CurrentResults {
 protected:
  const CompressedSubgraphs* input_;
  const uint32_t parent_index_;
  const TraverseOperator* owner_;
  ExpandEdgeTraverseContext* ctx_ = nullptr;

 public:
  CurrentResults(const CompressedSubgraphs* input, uint32_t parent_index, const TraverseOperator* owner,
                 ExpandEdgeTraverseContext* ctx)
      : input_(input), parent_index_(parent_index), owner_(owner), ctx_(ctx) {}

  virtual ~CurrentResults() {}

  virtual uint32_t getResults(uint32_t cap) = 0;
};

class ExpandEdgeSetToKeyTraverseContext : public ExpandEdgeTraverseContext {
  const CandidateSetView* candidates_ = nullptr;
  std::unique_ptr<CurrentResults> current_results_ = nullptr;
  uint64_t candidates_neighbor_size_ = 0;

 public:
  ExpandEdgeSetToKeyTraverseContext(const CandidateSetView* candidates, const void* graph,
                                    std::vector<CompressedSubgraphs>* outputs, QueryType profile,
                                    uint64_t candidate_neighbor_size = 0)
      : ExpandEdgeTraverseContext(candidates, graph, outputs, profile, false, nullptr),
        candidates_(candidates),
        candidates_neighbor_size_(candidate_neighbor_size) {}

  std::unique_ptr<TraverseContext> clone() const override {
    return std::make_unique<ExpandEdgeSetToKeyTraverseContext>(candidates_, getDataGraph<void>(), getOutputs(),
                                                               getQueryType(), candidates_neighbor_size_);
  }

  inline const CandidateSetView* getCandidateSet() const { return candidates_; }
  inline uint64_t getCandidateNeighborSize() const { return candidates_neighbor_size_; }
  inline bool hasRemainingResults() const { return current_results_ != nullptr; }
  inline uint32_t getResults(uint32_t cap) { return current_results_->getResults(cap); }
  inline void setResults(CurrentResults* ptr) { current_results_.reset(ptr); }

  template <typename G>
  bool init(const CandidateSetView& candidates, LabelID parent_label) {
    auto data_graph = getDataGraph<G>();
    candidates_neighbor_size_ = 0;
    for (auto candidate : candidates) {
      candidates_neighbor_size_ += data_graph->getVertexInDegreeWithHint(candidate, parent_label, 0);
    }
    if (candidates_neighbor_size_ == 0) {
      LOG(WARNING) << "WARNING: no neighbor for all " << candidates.size() << " candidates";
      return true;
    }
    return false;
  }
};

template <QueryType profile, typename G>
class CurrentResultsByCandidate : public CurrentResults {
 private:
  CandidateSetView::ConstIterator candidate_iter_;
  CandidateSetView::ConstIterator candidate_end_;
  const LabelID parent_label_;

 public:
  CurrentResultsByCandidate(const CompressedSubgraphs* input, uint32_t parent_index, LabelID parent_label,
                            const TraverseOperator* owner, ExpandEdgeSetToKeyTraverseContext* ctx)
      : CurrentResults(input, parent_index, owner, ctx),
        candidate_iter_(ctx->getCandidateSet()->begin()),
        candidate_end_(ctx->getCandidateSet()->end()),
        parent_label_(parent_label) {}

  uint32_t getResults(uint32_t cap) override {
    uint32_t n = 0;
    auto& parent_set = input_->getSet(parent_index_);

    auto exceptions = input_->getExceptions(owner_->getSameLabelKeyIndices(), owner_->getSameLabelSetIndices());
    auto data_graph = ctx_->getDataGraph<G>();
    for (; n < cap && candidate_iter_ != candidate_end_; ++candidate_iter_) {
      std::vector<VertexID> parents;
      VertexID candidate = *candidate_iter_;
      if (exceptions.count(candidate)) {
        continue;
      }
      if (owner_->filterTarget(candidate, *input_)) {
        continue;
      }
      auto neighbors = data_graph->getInNeighborsWithHint(candidate, parent_label_, 0);
      intersect(*parent_set, neighbors, &parents);  // No need for exceptions
      if
        constexpr(isProfileMode(profile)) {
          ctx_->updateIntersectInfo(parent_set->size() + neighbors.size(), parents.size());
        }
      if (parents.empty()) {
        continue;
      }
      auto output = ctx_->newOutput(*input_, parent_index_, newVertexSet(parents), candidate,
                                    owner_->getSameLabelSetIndices(), owner_->getSetPruningThreshold(), true);
      if (output == nullptr) continue;
#ifdef USE_FILTER
      if (owner_->filter(*output)) {
        ctx_->popOutput();
        continue;
      }
#endif
      ++n;
    }
    return n;
  }
};

template <QueryType profile, typename G, bool intersect_candidates>
class CurrentResultsByParent : public CurrentResults {
  unordered_set<VertexID> exceptions_;

 public:
  CurrentResultsByParent(const CompressedSubgraphs* input, uint32_t parent_index, const TraverseOperator* owner,
                         ExpandEdgeSetToKeyTraverseContext* ctx)
      : CurrentResults(input, parent_index, owner, ctx) {
    exceptions_ = input_->getExceptions(owner_->getSameLabelKeyIndices(), owner_->getSameLabelSetIndices());
  }

  uint32_t getResults(uint32_t cap) override {
    auto& parent_set = *(input_->getSet(parent_index_));
    unordered_map<VertexID, int> group_index;
    uint32_t old_count = ctx_->getOutputSize();
    auto target_label = owner_->getTargetLabel();
    auto& candidates = *((ExpandEdgeSetToKeyTraverseContext*)ctx_)->getCandidateSet();
    auto graph = ctx_->getDataGraph<G>();
    for (uint32_t i = 0; i < parent_set.size(); ++i) {
      auto parent_match = parent_set[i];
      if (exceptions_.count(parent_match)) continue;
      std::vector<VertexID> targets;
      auto neighbors = graph->getOutNeighborsWithHint(parent_match, target_label, 0);
      if (isProfileCandidateSIEffect(profile)) {
        intersectCandidateSetWithProfile(candidates, neighbors, &targets, exceptions_, ctx_, *input_,
                                         owner_->getTargetFilter());
      } else {
        if (intersect_candidates) {
          intersect(candidates, neighbors, &targets, exceptions_);
          if
            constexpr(isProfileMode(profile)) {
              ctx_->updateIntersectInfo(candidates.size() + neighbors.size(), targets.size());
            }
        } else {
          degreeFilter(neighbors, owner_->getTargetDegree(), graph, &targets, exceptions_);
        }
        if (!targets.empty()) {
          owner_->filterTargets(&targets, *input_);
        }
      }
      for (auto target : targets) {
        auto pos = group_index.find(target);
        if (pos == group_index.end()) {
          group_index[target] = ctx_->getOutputSize();
          auto output = ctx_->newOutput(*input_, parent_index_, singleVertexSet(parent_match), target,
                                        owner_->getSameLabelSetIndices(), owner_->getSetPruningThreshold(), false);
          if (output == nullptr) {
            group_index[target] = -1;
          }
        } else if (pos->second != -1) {  // the target is invalid if the group index is -1
          (*ctx_->getOutputs())[pos->second].UpdateSet(parent_index_, parent_match);
        }
      }
    }
#ifdef USE_FILTER
    auto n_pruned = owner_->filter((*ctx_->getOutputs()), old_count, ctx_->getOutputSize());
    ctx_->popOutputs(n_pruned);
#endif
    return ctx_->getOutputSize() - old_count;
  }
};  // class CurrentResultsByParent

template <QueryType profile, typename G, bool intersect_candidates>
class CurrentResultsByExtension : public CurrentResults {
 private:
  unordered_set<VertexID> seen_extensions_;
  std::queue<VertexID> extensions_;
  uint32_t parent_match_index_ = 0;
  unordered_set<VertexID> current_exceptions_;
  const LabelID parent_label_;

 public:
  CurrentResultsByExtension(const CompressedSubgraphs* input, uint32_t parent_index, LabelID parent_label,
                            const TraverseOperator* owner, ExpandEdgeTraverseContext* ctx)
      : CurrentResults(input, parent_index, owner, ctx), parent_label_(parent_label) {
    current_exceptions_ = input_->getExceptions(owner->getSameLabelKeyIndices(), owner->getSameLabelSetIndices());
  }

  uint32_t getResults(uint32_t cap) override {
    uint32_t n = 0;
    auto& parent_set = *input_->getSet(parent_index_);
    auto g = ctx_->getDataGraph<G>();
    while (true) {
      // check existing extensions first
      while (!extensions_.empty()) {
        auto candidate = extensions_.front();
        extensions_.pop();
        std::vector<VertexID> parents;  // valid parents for current candidate
        auto neighbors = g->getInNeighborsWithHint(candidate, parent_label_, 0);
        intersect(parent_set, neighbors, &parents);  // no need for exceptions
        if
          constexpr(isProfileMode(profile)) {
            ctx_->updateIntersectInfo(parent_set.size() + neighbors.size(), parents.size());
          }

        auto output = ctx_->newOutput(*input_, parent_index_, newVertexSet(parents), candidate,
                                      owner_->getSameLabelSetIndices(), owner_->getSetPruningThreshold(), true);
        if (output == nullptr) {
          continue;
        }
#ifdef USE_FILTER
        if (owner_->filter(*output)) {
          ctx_->popOutput();
          continue;
        }
#endif
        if (++n == cap) return n;
      }
      // all parent match extended
      if (parent_match_index_ == parent_set.size()) break;
      // get more extensions by extending from the next parent match
      auto parent_match = parent_set[parent_match_index_];
      getExtensions(parent_match, g);
      ++parent_match_index_;
    }
    return n;
  }

 private:
  inline void getExtensions(VertexID parent_match, const G* g) {
    if (current_exceptions_.count(parent_match)) {
      return;
    }
    std::vector<VertexID> current_extensions;
    auto neighbors = g->getOutNeighborsWithHint(parent_match, owner_->getTargetLabel(), 0);
    auto candidates = ((ExpandEdgeSetToKeyTraverseContext*)ctx_)->getCandidateSet();
    if (isProfileCandidateSIEffect(profile) && candidates != nullptr) {
      intersectCandidateSetWithProfile(*candidates, neighbors, &current_extensions, current_exceptions_, ctx_, *input_,
                                       owner_->getTargetFilter());
    } else {
      if (intersect_candidates) {
        intersect(*candidates, neighbors, &current_extensions, current_exceptions_);
        if
          constexpr(isProfileMode(profile)) {
            ctx_->updateIntersectInfo(candidates->size() + neighbors.size(), current_extensions.size());
          }
      } else {
        degreeFilter(neighbors, owner_->getTargetDegree(), g, &current_extensions, current_exceptions_);
      }
      if (!current_extensions.empty()) {
        owner_->filterTargets(&current_extensions, *input_);
      }
    }
    for (VertexID neighbor : current_extensions) {
      if (seen_extensions_.insert(neighbor).second) {
        extensions_.push(neighbor);
      }
    }
  }
};

template <typename G, bool intersect_candidates>
class ExpandEdgeSetToKeyOperator : public ExpandEdgeOperator {
  enum ExecutionMode { ByCandidate, ByParent, ByExtension };

 public:
  CONSTRUCT(SetToKey)

  uint32_t expand(uint32_t cap, TraverseContext* ctx) const override {
    return expandInner<QueryType::Execute>(cap, ctx);
  }

  uint32_t expandAndProfileInner(uint32_t cap, TraverseContext* ctx) const override {
    if (ctx->getQueryType() == QueryType::Profile) {
      return expandInner<QueryType::Profile>(cap, ctx);
    }
    if (ctx->getQueryType() == QueryType::ProfileCandidateSIEffect) {
      return expandInner<QueryType::ProfileCandidateSIEffect>(cap, ctx);
    }
    CHECK(ctx->getQueryType() == QueryType::ProfileWithMiniIntersection) << "unknown query type "
                                                                         << (uint32_t)ctx->getQueryType();
    return expandInner<QueryType::ProfileWithMiniIntersection>(cap, ctx);
  }

  std::string toString() const override {
    std::stringstream ss;
    ss << "ExpandEdgeSetToKeyOperator";
    toStringInner(ss);
    return ss.str();
  }

  std::pair<uint32_t, uint32_t> getOutputSize(const std::pair<uint32_t, uint32_t>& input_key_size) const override {
    return {input_key_size.first + 1, input_key_size.second + 1};
  }

  std::unique_ptr<TraverseContext> initTraverseContext(
      const CandidateSetView* candidates, std::vector<CompressedSubgraphs>* outputs, const void* graph,
      QueryType profile, const unordered_set<VertexID>* candidate_hashmap) const override {
    auto ret = std::make_unique<ExpandEdgeSetToKeyTraverseContext>(candidates, graph, outputs, profile);
    if (candidates != nullptr && ret->init<G>(*candidates, parent_label_)) {
      return nullptr;  // prune the whole traverse chain
    }
    return ret;
  }

 private:
  /* Assume the set intersection cost of two sorted sets of size n and m is 2(n+m).
   *
   * cost of enumerating candidates = total set intersection cost
   *                                = 2 * candidates_neighbor_size + 2 * |parent_set| * |candidates_|
   * cost of enumerating parent set
   *     = total set intersection cost + total CompressedSubgraphs lookup cost
   *     = 2 * set_neighbor_size_ + 2 * |parent_set| * |candidates_|
   *       + min(|parent_set| * |candidates_|, set_neighbor_size)
   */
  inline ExecutionMode getExecutionMode(const std::vector<VertexID>* parent_set, uint32_t cap,
                                        ExpandEdgeSetToKeyTraverseContext* ctx) const {
    if (ctx->getCandidateSet() == nullptr) return ByExtension;
    DCHECK_NE(ctx->getCandidateNeighborSize(), 0) << ctx->getCandidateSet()->size();
    uint64_t set_neighbor_size = 0;
    for (auto v : *parent_set) {
      set_neighbor_size += ctx->getDataGraph<G>()->getVertexOutDegreeWithHint(v, target_label_, 0);
    }
    auto enumerating_candidate_cost = 2 * ctx->getCandidateNeighborSize();
    auto enumerating_parent_cost =
        2 * set_neighbor_size + std::min(parent_set->size() * ctx->getCandidateSet()->size(), set_neighbor_size);
    if (enumerating_candidate_cost < enumerating_parent_cost) return ByCandidate;

    /* ByParent requires all parent to be processed before the outputs become ready, and thus the output size is unknown
     * in advance and unbounded. ByExtension incurs extra set intersection, but the output size can be bounded. */
    if (std::min(ctx->getCandidateSet()->size(), set_neighbor_size) <=
        cap) {             // here we use set_neighbor_size to approximate the output size. as long as the
                           // output size does not exceed cap, ByParent is preferred than ByExtension.
      return ByExtension;  // return ByParent;
    }
    return ByExtension;
  }

  template <QueryType profile>
  uint32_t expandInner(uint32_t cap, TraverseContext* base_ctx) const {
    auto ctx = dynamic_cast<ExpandEdgeSetToKeyTraverseContext*>(base_ctx);
    DCHECK_GT(cap, 0);
    uint32_t needed = cap;
    while (true) {
      if (ctx->hasRemainingResults()) {
        uint32_t got = ctx->getResults(needed);
        DCHECK_LE(got, needed);
        needed -= got;

        if (needed == 0) {
          return cap;
        }

        if
          constexpr(isProfileWithMiniIntersectionMode(profile)) {
            auto& parent_set = *ctx->getPreviousInput().getSet(parent_index_);
            for (auto parent_match : parent_set) {
              ctx->distinct_intersection_count += ctx->hasIntersectionParent(parent_match);
            }
          }
        ctx->setResults(nullptr);
      }
      if (!ctx->hasNextInput()) {
        return cap - needed;
      }
      expandInner<profile>(ctx->getCurrentInput(), needed, ctx);
      if
        constexpr(isProfileMode(profile)) {
          ctx->total_num_input_subgraphs += ctx->getCurrentInput().getNumSubgraphs();
        }
      ctx->nextInput();
    }
    return cap - needed;
  }

  template <QueryType profile>
  void expandInner(const CompressedSubgraphs& input, uint32_t cap, ExpandEdgeSetToKeyTraverseContext* ctx) const {
    auto& parent_set = input.getSet(parent_index_);
    ExecutionMode mode = getExecutionMode(parent_set.get(), cap, ctx);
    switch (mode) {
    case ByCandidate:
      // enumerating target candidates is likely to incur less computation
      ctx->setResults(new CurrentResultsByCandidate<profile, G>(&input, parent_index_, parent_label_, this, ctx));
      break;
    case ByParent:
      // enumerating vertices in the parent set is likely to incur less computation
      ctx->setResults(new CurrentResultsByParent<profile, G, intersect_candidates>(&input, parent_index_, this, ctx));
      break;
    default:
      // enumerating vertices in the parent set is likely to incur less computation
      ctx->setResults(new CurrentResultsByExtension<profile, G, intersect_candidates>(&input, parent_index_,
                                                                                      parent_label_, this, ctx));
    }
  }
};

TraverseOperator* ExpandEdgeOperator::newExpandEdgeKeyToSetOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex, const unordered_map<QueryVertexID, uint32_t>& indices,
    const std::vector<uint32_t>& same_label_key_indices, const std::vector<uint32_t>& same_label_set_indices,
    uint64_t set_pruning_threshold, std::unique_ptr<SubgraphFilter>&& filter, GraphType graph_type,
    bool intersect_candidates) {
  DCHECK_GT(indices.count(parent_vertex), 0);
  DCHECK_GT(indices.count(target_vertex), 0);
  return newTraverseOp<ExpandEdgeKeyToSetOperator>(
      graph_type, intersect_candidates, indices.at(parent_vertex), indices.at(target_vertex), parent_vertex,
      target_vertex, same_label_key_indices, same_label_set_indices, set_pruning_threshold, std::move(filter));
}

TraverseOperator* ExpandEdgeOperator::newExpandEdgeKeyToKeyOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex, const unordered_map<QueryVertexID, uint32_t>& indices,
    const std::vector<uint32_t>& same_label_key_indices, const std::vector<uint32_t>& same_label_set_indices,
    uint64_t set_pruning_threshold, std::unique_ptr<SubgraphFilter>&& filter, GraphType graph_type,
    bool intersect_candidates) {
  DCHECK_GT(indices.count(parent_vertex), 0);
  DCHECK_GT(indices.count(target_vertex), 0);
  return newTraverseOp<ExpandEdgeKeyToKeyOperator>(
      graph_type, intersect_candidates, indices.at(parent_vertex), indices.at(target_vertex), parent_vertex,
      target_vertex, same_label_key_indices, same_label_set_indices, set_pruning_threshold, std::move(filter));
}

// tricky case: expand, look up group for each match of target, and copy parent to set for each group; or
// consider an alternative: enumerate each candidate of target, and do set intersection between the target
// neighbors and parent sets
TraverseOperator* ExpandEdgeOperator::newExpandEdgeSetToKeyOperator(
    QueryVertexID parent_vertex, QueryVertexID target_vertex, const unordered_map<QueryVertexID, uint32_t>& indices,
    const std::vector<uint32_t>& same_label_key_indices, const std::vector<uint32_t>& same_label_set_indices,
    uint64_t set_pruning_threshold, std::unique_ptr<SubgraphFilter>&& filter, GraphType graph_type,
    bool intersect_candidates) {
  DCHECK_GT(indices.count(parent_vertex), 0);
  DCHECK_GT(indices.count(target_vertex), 0);
  return newTraverseOp<ExpandEdgeSetToKeyOperator>(
      graph_type, intersect_candidates, indices.at(parent_vertex), indices.at(target_vertex), parent_vertex,
      target_vertex, same_label_key_indices, same_label_set_indices, set_pruning_threshold, std::move(filter));
}

}  // namespace circinus

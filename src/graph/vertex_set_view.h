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

#include <iterator>
#include <utility>
#include <vector>

#include "glog/logging.h"

#include "algorithms/search.h"
#include "graph/types.h"

namespace circinus {

class VertexSetView {
  using RangeSize = uint32_t;
  size_t size_ = 0;

 protected:
  std::vector<std::pair<const VertexID*, RangeSize>> ranges_;

 public:
  class ConstIterator {
   private:
    const std::vector<std::pair<const VertexID*, RangeSize>>* ranges_;
    uint32_t range_idx_ = 0;
    RangeSize idx_in_range_ = 0;

   public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = int64_t;
    using value_type = VertexID;
    using pointer = const value_type*;
    using reference = const value_type&;

    friend class VertexSetView;

    friend std::ostream& operator<<(std::ostream& ss, const ConstIterator& iter) {
      return ss << iter.ranges_->size() << ", " << iter.range_idx_ << ", " << iter.idx_in_range_;
    }

    explicit ConstIterator(const std::vector<std::pair<const VertexID*, RangeSize>>& ranges, bool end = false)
        : ranges_(&ranges), range_idx_(ranges.size() * end) {}

    ConstIterator(const std::vector<std::pair<const VertexID*, RangeSize>>& ranges, uint32_t range_idx,
                  RangeSize idx_in_range)
        : ranges_(&ranges), range_idx_(range_idx), idx_in_range_(idx_in_range) {}

    inline reference operator*() const { return (*ranges_)[range_idx_].first[idx_in_range_]; }
    inline pointer operator->() const { return &this->operator*(); }

    /** prefix increment */
    inline ConstIterator& operator++() {
      bool next_range = (++idx_in_range_ == (*ranges_)[range_idx_].second);
      range_idx_ += next_range;
      idx_in_range_ *= (1 - next_range);
      return *this;
    }

    /** postfix increment */
    inline ConstIterator operator++(int) {
      auto ret = *this;
      ++(*this);
      return ret;
    }

    friend inline bool operator==(const ConstIterator& a, const ConstIterator& b) {
      DCHECK(a.ranges_ == b.ranges_) << "The iterators are from different views";
      return a.idx_in_range_ == b.idx_in_range_ && a.range_idx_ == b.range_idx_;
    }

    friend inline bool operator!=(const ConstIterator& a, const ConstIterator& b) {
      DCHECK(a.ranges_ == b.ranges_) << "The iterators are from different views";
      return a.idx_in_range_ != b.idx_in_range_ || a.range_idx_ != b.range_idx_;
    }

    /** search from this iterator to end */
    inline ConstIterator getLowerBound(const ConstIterator& end, value_type v) const {
      if (*this == end) {  // if this is end, return
        return *this;
      }

      uint32_t target_range_idx = range_idx_;
      // linear search for range as the number of ranges are supposed to be small
      while (target_range_idx < end.range_idx_ &&
             v > (*ranges_)[target_range_idx].first[(*ranges_)[target_range_idx].second - 1]) {
        ++target_range_idx;
      }
      if (target_range_idx == (*ranges_).size()) {  // end of view
        return end;
      }
      // find lower bound index in range
      const VertexID* first = (*ranges_)[target_range_idx].first + (target_range_idx == range_idx_) * idx_in_range_;
      const VertexID* last =
          (*ranges_)[target_range_idx].first +
          (target_range_idx == end.range_idx_ ? end.idx_in_range_ : (*ranges_)[target_range_idx].second);
      return ConstIterator((*ranges_), target_range_idx,
                           lowerBound(first, last, v) - (*ranges_)[target_range_idx].first);
    }
  };

  VertexSetView() = default;
  VertexSetView(const VertexID* start, const VertexID* end) { addRange(start, end); }
  explicit VertexSetView(const std::vector<VertexID>& vec) { addRange(vec.data(), vec.data() + vec.size()); }
  VertexSetView(ConstIterator& start, ConstIterator& end) {
    // no data
    if (start == end) return;
    auto& ranges = (*start.ranges_);
    // only one range
    auto range_begin = ranges[start.range_idx_].first;
    if (start.range_idx_ == end.range_idx_) {
      if (end.idx_in_range_ == start.idx_in_range_) {
        return;
      }
      addOffsetSize(range_begin + start.idx_in_range_, end.idx_in_range_ - start.idx_in_range_);
      return;
    }
    // maybe multiple ranges
    addOffsetSize(range_begin + start.idx_in_range_, ranges[start.range_idx_].second - start.idx_in_range_);
    for (uint32_t range_idx = start.range_idx_ + 1; range_idx < end.range_idx_; ++range_idx) {
      addOffsetSize(ranges[range_idx].first, ranges[range_idx].second);
    }
    if (end.range_idx_ < ranges.size() && end.idx_in_range_ > 0) {  // not the end of old range
      addOffsetSize(ranges[end.range_idx_].first, end.idx_in_range_);
    }
  }

  inline void addOffsetSize(const VertexID* start, RangeSize size) {
    DCHECK_NE(size, 0);
    size_ += size;
    ranges_.emplace_back(start, size);
  }

  inline void addRange(const VertexID* start, const VertexID* end) {
    DCHECK(start != nullptr);
    DCHECK(end != nullptr);
    if (start >= end) return;
    auto size = std::distance(start, end);
    size_ += size;
    if (!ranges_.empty() && ranges_.back().first + ranges_.back().second == start) {
      ranges_.back().second += size;
    } else {
      ranges_.emplace_back(start, size);
    }
  }

  inline size_t size() const { return size_; }
  inline bool empty() const { return size_ == 0; }
  inline const auto& getRanges() const { return ranges_; }
  inline void clear() {
    ranges_.clear();
    size_ = 0;
  }

  inline ConstIterator begin() const { return ConstIterator(ranges_); }
  inline ConstIterator end() const { return ConstIterator(ranges_, true); }
  inline ConstIterator iter(size_t idx) const {
    uint32_t range_idx = 0;
    RangeSize idx_in_range = 0;
    size_t offset = 0;
    for (; range_idx < ranges_.size(); ++range_idx) {
      if (offset + ranges_[range_idx].second > idx) {
        idx_in_range = idx - offset;
        break;
      }
      offset += ranges_[range_idx].second;
    }
    return ConstIterator(ranges_, range_idx, idx_in_range);
  }

  inline VertexID front() const {
    DCHECK_GT(ranges_.size(), 0);
    return *ranges_.front().first;
  }

  inline VertexID back() const {
    DCHECK_GT(ranges_.size(), 0);
    return ranges_.back().first[ranges_.back().second - 1];
  }
};

class SingleRangeVertexSetView : public VertexSetView {
 public:
  using ConstIterator = const VertexID*;
  SingleRangeVertexSetView() {}
  SingleRangeVertexSetView(const VertexID* start, size_t size) {
    if (start == nullptr) return;
    if (size == 0) return;
    addOffsetSize(start, size);
  }

  SingleRangeVertexSetView(const VertexID* start, const VertexID* end) {
    if (std::distance(start, end) == 0) return;
    addOffsetSize(start, std::distance(start, end));
  }

  inline VertexID operator[](uint32_t index) const {
    DCHECK(!ranges_.empty());
    DCHECK_LT(index, ranges_.front().second);
    return *(begin() + index);
  }

  inline ConstIterator begin() const { return ranges_.empty() ? nullptr : ranges_.front().first; }
  inline ConstIterator end() const {
    return ranges_.empty() ? nullptr : (ranges_.front().first + ranges_.front().second);
  }
};

}  // namespace circinus

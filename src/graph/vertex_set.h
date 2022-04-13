#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "graph/types.h"
#include "graph/vertex_set_view.h"

namespace circinus {

#define newSharedVSet std::make_shared<std::vector<VertexID>>

class VertexSet {
 private:
  std::shared_ptr<std::vector<VertexID>> data_ = nullptr;
  SingleRangeVertexSetView view_;

 public:
  VertexSet() {}

  explicit VertexSet(std::vector<VertexID>&& data) : data_(newSharedVSet(std::move(data))) {
    view_.addRange(data_->data(), data_->data() + data_->size());
  }

  explicit VertexSet(std::shared_ptr<std::vector<VertexID>>&& data) : data_(std::move(data)) {
    view_.addRange(data_->data(), data_->data() + data_->size());
  }

  explicit VertexSet(const SingleRangeVertexSetView& view) : view_(view) {}

  // update range
  VertexSet(const VertexSet& base, uint32_t offset, uint32_t size)
      : data_(base.data_), view_(base.view_.begin() + offset, size) {
    CHECK_LE(base.view_.begin() + offset + size, base.view_.end());
  }

  // std::vector<VertexID>* operator->() { return data_.get(); }
  const SingleRangeVertexSetView* operator->() const { return &view_; }
  const SingleRangeVertexSetView& operator*() const { return view_; }
  SingleRangeVertexSetView& operator*() { return view_; }

  std::vector<VertexID>* get() const { return data_.get(); }

  const auto& getBuffer() const { return data_; }

  void push_back(VertexID v) {
    data_->push_back(v);
    view_ = SingleRangeVertexSetView(data_->data(), data_->size());
  }

  inline void resetTargetView() {
    if (data_->empty()) return;
    view_ = SingleRangeVertexSetView(data_->data(), data_->size());
  }

  inline std::vector<VertexID>& resetTargets() {
    view_.clear();
    data_ = newSharedVSet();
    return *data_;
  }

  inline void setTargetView(const std::shared_ptr<std::vector<VertexID>>& buffer,
                            const SingleRangeVertexSetView& view) {
    data_ = buffer;
    view_ = view;
  }
};

#define singleVertexSet(vertex) VertexSet(std::vector<VertexID>({vertex}))
#define newVertexSet(set) VertexSet(std::move(set))

}  // namespace circinus

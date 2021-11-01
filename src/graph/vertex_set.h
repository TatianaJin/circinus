// Copyright 2021 HDL
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
// the License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

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

  void push_back(VertexID v) {
    data_->push_back(v);
    view_ = SingleRangeVertexSetView(data_->data(), data_->size());
  }
};

#define singleVertexSet(vertex) VertexSet(std::vector<VertexID>({vertex}))
#define newVertexSet(set) VertexSet(std::move(set))

}  // namespace circinus

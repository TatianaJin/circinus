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
#include <vector>

#include "glog/logging.h"

#include "graph/types.h"
#include "graph/vertex_set_view.h"

namespace circinus {

class CandidateSetView : public VertexSetView {
 public:
  CandidateSetView(const std::vector<VertexID>* candidates, const CandidateScope& scope,
                   const std::vector<VertexID>& partition_offsets);
};

}  // namespace circinus

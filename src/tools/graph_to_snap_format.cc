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

#include <iostream>

#include "glog/logging.h"

#include "graph/types.h"
#include "utils/file_utils.h"

using circinus::EdgeID;
using circinus::LabelID;
using circinus::VertexID;

int main(int argc, char** argv) {
  if (argc == 1) {
    std::cout << "./GraphToSNAPFormat input [output]" << std::endl;
    return 0;
  }

  std::string input = argv[1];
  std::string output = (argc == 3) ? argv[2] : (input + ".snap");

  auto infile = circinus::openFile(input);
  auto outfile = circinus::openOutputFile(output);

  char line_type;
  VertexID n_vertices;
  EdgeID n_edges;

  // process line: t n_vertices n_edges
  infile >> line_type >> n_vertices >> n_edges;
  outfile << "# Nodes: " << n_vertices << " Edges: " << n_edges << std::endl;

  // the first n_vertices lines should be of type v, with continuous vertex id from 0
  VertexID id;
  LabelID label;
  VertexID degree;
  for (uint32_t i = 0; i < n_vertices; ++i) {
    CHECK(infile >> line_type >> id >> label >> degree);
  }

  // next n_edges lines should be of type e, output as edge list
  outfile << "# FromNodeId\tToNodeId" << std::endl;
  VertexID v1, v2;
  for (EdgeID i = 0; i < n_edges; ++i) {
    CHECK(infile >> line_type) << i + n_vertices + 2 << " " << input;
    if (line_type != 'e') {
      CHECK(infile >> line_type >> v1 >> v2) << i + n_vertices + 2 << " " << input;
    } else {
      CHECK(infile >> v1 >> v2) << i + n_vertices + 2 << " " << input;
    }
    CHECK_EQ(line_type, 'e');
    outfile << v1 << '\t' << v2 << std::endl;
  }
  infile.close();
  outfile.close();
}

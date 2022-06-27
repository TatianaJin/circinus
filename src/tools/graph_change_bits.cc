#include "glog/logging.h"

#include "graph/graph.h"
#include "graph/partitioned_graph.h"

using circinus::Graph;
using circinus::ReorderedPartitionedGraph;
using circinus::openFile;
using circinus::openOutputFile;
using circinus::VertexID;
using circinus::EdgeID;
using circinus::LabelID;
using circinus::unordered_map;

#define V_input_t uint64_t
#define E_input_t uint64_t

template <typename T>
static void vectorToBinaryStream(std::ostream& output, const std::vector<T>& vec) {
  size_t size = vec.size();
  output.write(reinterpret_cast<const char*>(&size), sizeof(size));
  output.write(reinterpret_cast<const char*>(vec.data()), sizeof(T) * vec.size());
}

template <typename T>
static void binaryStreamToVector(std::istream& input, std::vector<T>& vec) {
  size_t size;
  input.read(reinterpret_cast<char*>(&size), sizeof(size));
  vec.resize(size);
  input.read(reinterpret_cast<char*>(vec.data()), sizeof(T) * vec.size());
}

template <typename InputT, typename OutputT>
void inputToOutput(std::istream& in, std::ostream& out, std::vector<InputT>& in_buf, std::vector<OutputT>& out_buf) {
  size_t size;
  in.read(reinterpret_cast<char*>(&size), sizeof(size));
  out.write(reinterpret_cast<const char*>(&size), sizeof(size));
  in_buf.resize(size);
  in.read(reinterpret_cast<char*>(in_buf.data()), sizeof(InputT) * in_buf.size());
  out_buf.assign(in_buf.begin(), in_buf.end());
  out.write(reinterpret_cast<const char*>(out_buf.data()), sizeof(OutputT) * out_buf.size());
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  google::SetUsageMessage("./GraphChangeBits input_path [output_path]");
  if (FLAGS_log_dir == "") {
    google::LogToStderr();
  }

  if (argc == 1) {
    LOG(INFO) << gflags::ProgramUsage();
    return 0;
  }

  std::string input_file = argv[1];
  std::string output_file;
  if (argc == 3) {
    output_file = argv[2];
  } else {
    output_file = input_file + "." + std::to_string(sizeof(VertexID));
  }
  auto input = openFile(input_file, std::ios::binary);
  auto output = openOutputFile(output_file, std::ios::binary);

  // graph type flag
  bool partitioned_graph = false;
  input.read(reinterpret_cast<char*>(&partitioned_graph), sizeof(bool));
  output.write(reinterpret_cast<char*>(&partitioned_graph), sizeof(bool));
  LOG(INFO) << "partitioned_graph = " << partitioned_graph;

  /* input >>> */
  V_input_t n_vertices_ = 0;
  E_input_t n_edges_ = 0;
  V_input_t max_degree_ = 0;
  std::vector<E_input_t> vlist_;
  std::vector<V_input_t> elist_;
  /* <<< input */
  /* output >>> */
  VertexID n_vertices = 0;
  EdgeID n_edges = 0;
  VertexID max_degree = 0;
  std::vector<EdgeID> vlist;
  std::vector<VertexID> elist;
  /* <<< output */
  std::vector<char> buf;

  auto floaddump = [&]() {
    input.read(reinterpret_cast<char*>(&n_vertices_), sizeof(n_vertices_));
    input.read(reinterpret_cast<char*>(&n_edges_), sizeof(n_edges_));
    input.read(reinterpret_cast<char*>(&max_degree_), sizeof(max_degree_));
    LOG(INFO) << "Vertices " << n_vertices_ << ", edges " << n_edges_ << ", max degree " << max_degree_;
    n_vertices = n_vertices_;
    n_edges = n_edges_;
    max_degree = max_degree_;
    output.write(reinterpret_cast<char*>(&n_vertices), sizeof(n_vertices));
    output.write(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));
    output.write(reinterpret_cast<char*>(&max_degree), sizeof(max_degree));

    inputToOutput(input, output, vlist_, vlist);
    CHECK_EQ(vlist_.size(), n_vertices_ + 1);
    inputToOutput(input, output, elist_, elist);
    CHECK_EQ(elist_.size(), 2 * n_edges_);

    {  // vertex_cardinality_by_label_
      size_t map_size;
      input.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
      output.write(reinterpret_cast<char*>(&map_size), sizeof(map_size));
      buf.resize(map_size * (sizeof(LabelID) + sizeof(uint32_t)));
      input.read(buf.data(), buf.size());
      output.write(buf.data(), buf.size());
    }
  };

  size_t size;
  if (partitioned_graph) {
    uint32_t n_partitions;
    input.read(reinterpret_cast<char*>(&n_partitions), sizeof(n_partitions));
    output.write(reinterpret_cast<char*>(&n_partitions), sizeof(n_partitions));
    floaddump();
    uint64_t n_edge_cuts_ = 0;
    input.read(reinterpret_cast<char*>(&n_edge_cuts_), sizeof(n_edge_cuts_));
    output.write(reinterpret_cast<char*>(&n_edge_cuts_), sizeof(n_edge_cuts_));
    LOG(INFO) << "n_edge_cuts_ = " << n_edge_cuts_;

    inputToOutput(input, output, elist_, elist);  // vertex_ids_
    CHECK_EQ(elist_.size(), n_vertices_);
    inputToOutput(input, output, elist_, elist);  // partition_offsets_
    CHECK_EQ(elist_.size(), n_partitions + 1);

    {  // label_idx_
      input.read(reinterpret_cast<char*>(&size), sizeof(size));
      output.write(reinterpret_cast<char*>(&size), sizeof(size));
      buf.resize(size * sizeof(LabelID) * 2);
      input.read(buf.data(), buf.size());
      output.write(buf.data(), buf.size());
    }
    {  // label_set_
      input.read(reinterpret_cast<char*>(&size), sizeof(size));
      output.write(reinterpret_cast<char*>(&size), sizeof(size));
      buf.resize(size * sizeof(LabelID));
      input.read(buf.data(), buf.size());
      output.write(buf.data(), buf.size());
    }
    {  // label ranges per partition
      input.read(reinterpret_cast<char*>(&size), sizeof(size));
      output.write(reinterpret_cast<char*>(&size), sizeof(size));
      for (size_t j = 0; j < size; ++j) {
        inputToOutput(input, output, elist_, elist);
      }
    }
  } else {  // if (partitioned_graph)
    floaddump();
    input.read(reinterpret_cast<char*>(&size), sizeof(size));
    output.write(reinterpret_cast<char*>(&size), sizeof(size));
    buf.resize(size * sizeof(LabelID));
    input.read(buf.data(), buf.size());
    output.write(buf.data(), buf.size());
  }

  input.close();
  output.close();
  return 0;
}

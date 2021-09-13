#! /usr/bin/env python3

import atexit
import numpy as np
import os.path as osp
import pickle
import socket
import struct
import zmq

from argparse import ArgumentParser
from datetime import datetime
from networkx import Graph


def unpack(data, mode='str', num=1):
  if mode == "bool":
    return struct.unpack('?' * num, data)
  if mode == "double":
    return struct.unpack('d' * num, data)
  if mode == "uint64_t":
    return struct.unpack('Q' * num, data)
  if mode == "uint32_t":
    return struct.unpack('I' * num, data)
  if mode == "str":
    return (struct.unpack('{0}s'.format(len(data)), data)[0].decode('utf-8'))


def parse_args():
  now = datetime.now().strftime("%Y%m%d%H%M%S")
  parser = ArgumentParser("Circinus Cost Learner")
  group = parser.add_mutually_exclusive_group()
  group.add_argument("-d",
                     "--data_path",
                     default=osp.join(osp.dirname(osp.dirname(osp.realpath(__file__))), "cost_learner_data_" + now),
                     help="The path to store historical records")
  group.add_argument("--history", help="The path to retrieve historical records")
  return parser.parse_args()


class QueryRecord:

  def __init__(self, n_labels, n_partitions, n_qvs, pqv, vertex_labels=None, nx_graph=None):
    self.n_labels = n_labels
    self.n_partitions = n_partitions
    self.n_qvs = n_qvs
    self.pqv = pqv
    self.vertex_labels = np.array(vertex_labels)
    self.graph = nx_graph
    self.plans = []
    self.n_records = 0

  def add_plan(self, plan_features, record_features, record_truths, record_remove_edges=None):
    assert (len(record_features) == len(record_truths)), "record features {0} truths {1}".format(len(record_features), len(record_truths))
    assert (record_features[0].shape[0] == plan_features.shape[0]), "n_vertex in record {0}, in plan {1}".format(
      record_features[0].shape[0], plan_features.shape[0])
    self.plans.append((plan_features, record_features, record_truths, record_remove_edges))
    self.n_records += len(record_truths)

  def get_records(self):
    query_features = np.eye(self.n_labels, dtype=int)[self.vertex_labels]
    for plan_features, record_features, record_truths, record_remove_edges in self.plans:
      for record_i in range(len(record_truths)):
        h = np.concatenate((query_features, plan_features, record_features[record_i]), axis=1)
        graph = self.graph
        if record_remove_edges is not None and record_i in record_remove_edges:
          graph = self.graph.copy()
          for edges in record_remove_edges[record_i]:
            graph.remove_edge(*edges)
        yield graph, h, record_truths[record_i]


class CircinusCostLearner:

  def __init__(self, args):
    self.context = zmq.Context()
    self.update_listen_sock = zmq.Socket(self.context, zmq.PULL)
    self.update_listen_address = "tcp://{0}:{1}".format(socket.gethostname(), self.update_listen_sock.bind_to_random_port("tcp://*"))
    print("Update listener on {0}".format(self.update_listen_address))
    self.record_store_path = args.data_path
    self.query_records = []
    atexit.register(self.save_records)

  def save_records(self):
    if len(self.query_records) == 0:
      return
    print("Saving records to", self.record_store_path)
    with open(self.record_store_path, 'wb') as f:
      pickle.dump(self.query_records, f)
    print("Records saved at", self.record_store_path)

  def _listen_update(self):
    msgs = self.update_listen_sock.recv_multipart()
    n_labels = unpack(msgs[0], 'uint32_t')[0]
    n_partitions = unpack(msgs[1], 'uint32_t')[0]
    n_qvs = unpack(msgs[2], 'uint32_t')[0]
    vertex_labels = unpack(msgs[3], 'uint32_t', n_qvs)
    query = Graph()
    for v in range(n_qvs):
      query.add_node(v)
    msg_i = 4
    for v in range(n_qvs):
      msg_i, n_v_nbrs = msg_i + 1, unpack(msgs[msg_i], 'uint32_t')[0]
      msg_i, adj = msg_i + 1, unpack(msgs[msg_i], 'uint32_t', n_v_nbrs)
      for nb in adj:
        if nb > v:
          query.add_edge(v, nb)
    msg_i, pqv = msg_i + 1, unpack(msgs[msg_i], 'uint32_t')[0]
    qr = QueryRecord(n_labels, n_partitions, n_qvs, pqv, vertex_labels=vertex_labels, nx_graph=query)

    msg_i, n_logical_plans = msg_i + 1, unpack(msgs[msg_i], 'uint32_t')[0]
    for plan_i in range(n_logical_plans):
      msg_i, vertex_cardinality = msg_i + 1, unpack(msgs[msg_i], 'uint32_t', n_qvs)
      msg_i, vertex_partitions = msg_i + 1, unpack(msgs[msg_i], 'str')[0]
      msg_i, matching_order = msg_i + 1, unpack(msgs[msg_i], 'uint32_t', n_qvs)
      msg_i, n_records = msg_i + 1, unpack(msgs[msg_i], 'uint32_t')[0]
      vertex_partition_features = np.ones((n_qvs, n_partitions), dtype=int)
      vertex_partition_features[pqv] = np.zeros_like(vertex_partition_features[pqv], dtype=int)
      for part in vertex_partitions.split(','):
        if part != "":
          vertex_partition_features[pqv][int(part)] = 1
      record_features = [None] * n_records
      record_remove_edges = None
      for record_i in range(n_records):
        msg_i, subquery_size = msg_i + 1, unpack(msgs[msg_i], 'uint32_t')[0]
        msg_i, cover_bits = msg_i + 1, unpack(msgs[msg_i], 'uint64_t')[0]
        msg_i, remove_parent_edge_flag = msg_i + 1, unpack(msgs[msg_i], 'bool')[0]
        subquery_mask_feature = np.zeros((n_qvs, 2), dtype=int)
        existing_vertices = set()
        for i in range(subquery_size):
          subquery_mask_feature[matching_order[i]][0] = 1
          subquery_mask_feature[matching_order[i]][1] = (cover_bits >> matching_order[i]) & 1
          existing_vertices.add(matching_order[i])
        record_features[record_i] = subquery_mask_feature
        if remove_parent_edge_flag:  # remove set-parent to target edges
          remove_edges = []
          for nb in query.adj[matching_order[subquery_size - 1]]:
            if nb in existing_vertices and (cover_bits >> nb) & 1 == 0:
              remove_edges.append((matching_order[subquery_size - 1], nb))
          if record_remove_edges is None:
            record_remove_edges = {record_i: remove_edges}
          else:
            record_remove_edges.update({record_i: remove_edges})
      record_truths = [None] * n_records
      for record_i in range(n_records):
        msg_i, truth = msg_i + 1, unpack(msgs[msg_i], 'uint64_t')[0]
        record_truths[record_i] = float(truth)
      plan_features = np.concatenate((np.array(vertex_cardinality).reshape(n_qvs, 1), vertex_partition_features), axis=1)
      qr.add_plan(plan_features, record_features, record_truths, record_remove_edges=record_remove_edges)
    self.query_records.append(qr)

  def get_records(self):
    for record in self.query_records:
      for entry in record.get_records():
        yield entry

  def serve(self):
    while True:
      self._listen_update()


def load_records(path):
  with open(path, 'rb') as f:
    return pickle.load(f)


if __name__ == '__main__':
  args = parse_args()
  if args.history:
    query_records = load_records(args.history)
    for record in query_records:
      for entry in record.get_records():
        print(entry)
  else:
    learner = CircinusCostLearner(args)
    learner.serve()

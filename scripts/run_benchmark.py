#! /usr/bin/env python3

import argparse
import os
from os import path as osp
from sys import stderr


def get_args():
  parser = argparse.ArgumentParser(description='Benchmark')
  parser.add_argument('-p', '--project_dir', help='Circinus project dir.', required=True)
  parser.add_argument('-t', '--time_out', default=300, help='Query timeout.')
  parser.add_argument('-m', '--match', default='Benchmark', help='The executable for running circinus.')

  # workload
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('-c', '--config', help='The file of query configs (under project_dir/evaluation/configs).')
  group.add_argument('-q', '--query', help='A single line of query config.')

  # opts for workload
  parser.add_argument('-l', '--limit', default=0, help='Match limit. 0 means unlimited.')
  parser.add_argument('-b', '--batch', default=1024, help='Batch size.')
  parser.add_argument('-f', '--filter', default='cfl', help='Candidate filter.')
  parser.add_argument('-o', '--order', default='cfl', help='Matching order strategy.')
  parser.add_argument('-s', '--strategy', default='dynamic', choices=['dynamic', 'eager', 'static', 'all'], help='Vertex cover strategy.')
  parser.add_argument('--profile', choices=['simple', 'si'])
  parser.add_argument('--partition', type=int, default=1, help='The number of partitions for graph')
  parser.add_argument('--upg', '--use_partitioned_graph', type=int, default=1, choices=[0, 1])
  parser.add_argument('--ipp', '--intra_partition_plan', type=int, default=1, choices=[0, 1])
  parser.add_argument('--pqv', default='none', choices=['none', 'cc'])

  args = parser.parse_args()
  args.bin_dir = osp.join(args.project_dir, "build", "tests")
  args.log_dir = osp.join(args.project_dir, "evaluation", "log")
  args.config_dir = osp.join(args.project_dir, "evaluation", "configs")
  return args, parser


def get_log_path(args):
  base = "_".join([
    str(x) for x in (args.match, args.limit, args.batch, args.filter, args.order, args.strategy, "partition{0}".format(args.partition),
                     "upg{0}".format(args.upg), "ipp{0}".format(args.ipp), args.pqv, args.config) if x is not None
  ])
  path = osp.join(args.log_dir, base)
  if args.profile is not None:  # for profiling, we output the profile for each query separately in a folder with the result log
    log_dir = "{0}_{1}".format(path, args.profile)
    if not osp.exists(log_dir):
      os.mkdir(log_dir)
    log = osp.join(log_dir, "log")
    return log_dir, log
  return None, path


def get_common_flags(args, log):
  return "-match_limit {0} -batch_size {1} -filter {2} -vertex_cover {3} -output_file {4} -match_order {5} -upg {upg} -ipp {ipp} -pqv {pqv} -partition {partition}".format(
    args.limit, args.batch, args.filter, args.strategy, log, args.order, upg=args.upg, ipp=args.ipp, pqv=args.pqv, partition=args.partition)


def run_batch(args):
  config_path = osp.join(args.config_dir, args.config)
  if not osp.isfile(config_path):
    stderr.write("File not exists {0}\n".format(config_path))
    return

  log_dir, log = get_log_path(args)

  with open(log, 'a') as log_f:
    log_f.write("dataset,query_size,query_mode,query_index,elapsed_time,filter_time,plan_time,enumerate_time,n_embeddings,order,max_task_time\n")

  common_flags = get_common_flags(args, log)
  with open(config_path, 'r') as config_f:
    for line in config_f:
      if line.startswith("dataset") or line.startswith("#"):
        continue
      run_query(line.strip(), args, log_dir, common_flags)


def run_query(config, args, log_dir, common_flags):
  flags = config.split(',')
  profile_flag = "" if args.profile is None else "-profile_prefix {1} -profile {0}".format(1 if args.profile == "simple" else 2, log_dir)
  cmd = "timeout {0} {1} -verbosity 0 -dataset {2} -query_size {3} -query_mode {4} -query_index {5} {6} {7}".format(
    args.time_out, osp.join(args.bin_dir, args.match), flags[0], flags[1], flags[2], flags[3], profile_flag, common_flags)
  os.system(cmd)


if __name__ == '__main__':
  args, _ = get_args()
  if not osp.exists(args.log_dir):
    os.mkdir(args.log_dir)
  print(args)
  if args.config is not None:
    run_batch(args)
  else:
    log_dir, log = get_log_path(args)
    run_query(args.query, args, log_dir, get_common_flags(args, log))

#! /usr/bin/env python3

import argparse
import pandas as pd

from os import path as osp
from sys import stderr


def get_args():
  parser = argparse.ArgumentParser(description='Experimental analysis')
  parser.add_argument('target_file', help='The file of results to analyze')
  parser.add_argument('--check', help='The groundtruth file for checking whether the n_embeddings values in the target file are correct')
  parser.add_argument('-t', '--time', help='The baseline file for comparing enumerate time')
  parser.add_argument('-m', '--missing', help='The config file for getting missing queries')
  parser.add_argument('--config', help='Config output path')

  args = parser.parse_args()
  return args, parser


def read(path):
  f = pd.read_csv(path)
  return f.set_index(['dataset', 'query_size', 'query_mode', 'query_index'])


def check_embeddings(target, groundtruth_file):
  """ Returns the number of wrong n_embeddings results """
  groundtruth = read(groundtruth_file)[["n_embeddings"]].dropna()
  comparison = groundtruth.join(target[["n_embeddings"]], lsuffix='_truth', rsuffix='_check').dropna()
  mismatch = comparison.loc[comparison['n_embeddings_truth'] != comparison['n_embeddings_check']]
  if len(mismatch) > 0:
    print(mismatch)
    return len(mismatch)
  return 0


def compare_enumerate_time(target, baseline_file):
  baseline = read(baseline_file)[["enumerate_time"]].dropna()
  print("Complete queries (target/baseline):\t{0}/{1}".format(len(target), len(baseline)))
  comparison = baseline.join(target[['enumerate_time']], lsuffix='_base', rsuffix='_check').dropna()
  speedup = comparison['enumerate_time_base'] - comparison['enumerate_time_check']
  comparison['speedup'] = speedup
  comparison['speedup_ratio'] = speedup / comparison['enumerate_time_base']
  print(comparison[['speedup', 'speedup_ratio']].describe())
  if sum(speedup < 0) > 0:
    print(comparison[["speedup", 'enumerate_time_base', 'speedup_ratio']].loc[speedup < 0])


def get_missing_query_config(target, config_file):
  baseline = read(config_file)[["enumerate_time"]]
  comparison = baseline.join(target[['enumerate_time']], lsuffix='', rsuffix='_check')
  missing = comparison.loc[comparison['enumerate_time_check'].isnull()][['enumerate_time']]
  missing.to_csv("{0}_missing.csv".format(config_file))
  print("Missing queries are output to {0}_missing.csv".format(config_file))


def get_config(target, config_file):
  target = target[[]]
  target.to_csv(config_file)
  print("Query config is output to {0}".format(config_file))


if __name__ == '__main__':
  args, _ = get_args()
  target = read(args.target_file).dropna()
  if args.check:
    print("Mismatch:\t{0}".format(check_embeddings(target, args.check)))
  if args.time:
    compare_enumerate_time(target, args.time)
  if args.missing:
    print("Total elapsed backtracking time", target['elapsed_time'].sum())
    get_missing_query_config(target, args.missing)
  if args.config:
    get_config(target, args.config)

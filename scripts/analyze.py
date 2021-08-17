#! /usr/bin/env python3

import argparse
import pandas as pd

from os import path as osp
from os import listdir
from sys import stderr


def get_args():
  parser = argparse.ArgumentParser(description='Experimental analysis')
  parser.add_argument('target_file', help='The file of results to analyze, or the profile directory')
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
  print(speedup.sum())
  print(comparison.iloc[speedup.argmin()])
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


def read_profile(profile_path, df):
  profile = []
  config = osp.basename(profile_path)
  with open(profile_path, 'r') as f:
    n_ops = None
    step = 1
    for line in f:
      line = line.strip()
      if line.startswith("[ Plan"):
        n_ops = int(line.split(" ")[-2])
        #print("#ops = {0}".format(line.split(" ")[-2]))
      elif n_ops is not None:
        splits = line.split(',')
        if int(splits[0]) == n_ops - 1:
          n_ops = None
        elif int(splits[0]) != 0:
          op = splits[1].split(' ', 2)[0]
          if op == "ExpandIntoOperator":
            profile[-1][2] += '-Into'
            profile[-1][4] += int(splits[7]) # si_count
            profile[-1][4] += int(splits[8]) # si_input
            profile[-1][5] += float(splits[2]) # time
          else:
            profile.append([config, step, op, 0.0, float(splits[7]), float(splits[8]), float(splits[2])])
            step = step + 1
      elif line.startswith("step_costs"):
        for i, v in enumerate(line.split(' ')[2:]):
            profile[i][3] = float(v)
  update = pd.DataFrame(data=profile, columns=['query', 'step', 'op', 'cost', 'si_count', 'si_input', 'time'])
  #print(update)
  return df.append(update)


def analyze_profile(folder):
  df = pd.DataFrame(columns=['query', 'step', 'op', 'cost', 'si_count', 'si_input', 'time'])
  for name in listdir(folder):
    if name == "log" or name.endswith(".swp"):
      continue
    df = read_profile(osp.join(folder, name), df)
  print(df)
  print("================================== Correlation ==================================")
  print(df[['cost', 'si_count', 'si_input', 'time']].corr())
  print("=============================== Correlation By Op ===============================")
  print(df.groupby('op').corr(min_periods=1))
  print("============================== Correlation By Step ==============================")
  print(df.groupby('step').corr(min_periods=1))


if __name__ == '__main__':
  args, _ = get_args()
  if osp.isdir(args.target_file):
    target_log = osp.join(args.target_file, "log")
    analyze_profile(args.target_file)
  else:
    target_log = args.target_file
  target = read(target_log).dropna()
  if args.check:
    print("Mismatch:\t{0}".format(check_embeddings(target, args.check)))
  if args.time:
    compare_enumerate_time(target, args.time)
  if args.missing:
    print("Total elapsed backtracking time", target['elapsed_time'].sum())
    get_missing_query_config(target, args.missing)
  if args.config:
    get_config(target, args.config)

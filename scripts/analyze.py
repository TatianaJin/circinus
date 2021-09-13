#! /usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

from os import path as osp
from os import listdir
from sklearn.linear_model import LinearRegression, ElasticNet, Ridge
from sklearn.model_selection import KFold
from sklearn.neural_network import MLPRegressor
from scipy.sparse import csr_matrix
from sys import stderr
from tqdm import tqdm


def get_args():
  parser = argparse.ArgumentParser(description='Experimental analysis')
  parser.add_argument('target_file', help='The file of results to analyze, or the profile directory')
  parser.add_argument('--check', help='The groundtruth file for checking whether the n_embeddings values in the target file are correct')
  parser.add_argument('-t', '--time', help='The baseline file for comparing enumerate time')
  parser.add_argument('-m', '--missing', help='The config file for getting missing queries')
  parser.add_argument('--config', help='Config output path')
  parser.add_argument('--parallel', action='store_true', help='Parallel execution analysis')
  parser.add_argument('--plan', help='The baseline file for comparing plan time')
  parser.add_argument('--topo', help='The file path of query graph topological features')

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


def compare_time(target, baseline_file, time_name=['enumerate_time']):
  if time_name is None or len(time_name) == 0:
      return
  baseline = read(baseline_file)[time_name].dropna()
  if len(time_name) == 1:
      print("Complete queries (target/baseline):\t{0}/{1}".format(len(target), len(baseline)))
      target = target[time_name]
      comparison_name = time_name[0]
  else:
      print("Complete queries (target/baseline):\t{0}/{1}".format(len(target), len(baseline)))
      baseline['sum_time'] = baseline[time_name[0]]
      for col in time_name[1:]:
          baseline['sum_time'] += baseline[col]
      #baseline = baseline[['sum_time']]
      target['sum_time'] = target[time_name[0]]
      for col in time_name[1:]:
          target['sum_time'] += target[col]
      #target = target[['sum_time']]
      comparison_name = 'sum_time'
  comparison = baseline.join(target, lsuffix='_base', rsuffix='_check').dropna()
  speedup = comparison[comparison_name + '_base'] - comparison[comparison_name + '_check']
  comparison['speedup'] = speedup
  comparison['speedup_ratio'] = speedup / comparison[comparison_name + '_base']
  print(comparison[['speedup', 'speedup_ratio']].describe())
  print("Total speedup (seconds)", speedup.sum())
  # print(comparison.iloc[speedup.argmin()])
  if sum(speedup < 0) > 0:
    print(comparison[["speedup", comparison_name + '_base', 'speedup_ratio']].loc[speedup < 0].sort_values(by='speedup', ascending=True))


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


def parallel_analysis(target):
  #elapsed = 'elapsed_execution_time'
  elapsed = 'elapsed_time'
  target = target[['enumerate_time', 'max_task_time', elapsed]]
  print("========== parallel ratio = enumerate_time / elapsed_time ==========")
  print((target['enumerate_time'] / target[elapsed]).describe())
  print("======= dominant task ratio = max_task_time / enumerate_time =======")
  max_task_ratio = target['max_task_time'] / target['enumerate_time']
  print(max_task_ratio.describe())
  max_task_ratio.nlargest(10)


def read_profile(profile_path, df):
  profile = []
  config = osp.basename(profile_path)
  with open(profile_path, 'r') as f:
    n_ops = None
    for line in f:
      line = line.strip()
      if line.startswith("[ Plan"):
        step = 1
        n_ops = int(line.split(" ")[-2])
        #print("#ops = {0}".format(line.split(" ")[-2]))
      elif n_ops is not None:
        splits = line.split(',')
        if int(splits[0]) == n_ops - 1:
          n_ops = None
        elif int(splits[0]) != 0:
          op, detail = splits[1].split(' ', 1)
          parents, _ = detail.split(' -> ')
          n_parents = len(parents.split(' '))
          if op == "ExpandIntoOperator":
            profile[-1][2] += '-Into'
            profile[-1][4] += int(splits[7])  # si_count
            profile[-1][5] += int(splits[8])  # si_input
            profile[-1][6] += float(splits[2])  # time
            profile[-1][8] += "-" + str(n_parents)
          else:
            profile.append([
              config, step, op, 0.0,
              float(splits[7]),
              float(splits[8]),
              float(splits[2]), 'dense' if 'dense' in config else 'sparse',
              str(n_parents),
              int(splits[-1])
            ])
            step = step + 1
      elif line.startswith("step_costs"):
        costs = line.split(' ')[2:]
        offset = len(profile) - len(costs)
        for i, v in enumerate(costs):
          profile[i + offset][3] = float(v)
  update = pd.DataFrame(data=profile,
                        columns=['query', 'step', 'op', 'cost', 'si_count', 'si_input', 'time', 'mode', 'n_parents', 'candidate_si_effect'])
  #print(update)
  #exit(0)
  return df.append(update) if df is not None else update


def topo_performance_analysis(target, topo_path, baseline_file):
  query_configs = []
  label_sequence_map = dict()
  indptr = [0]  # csr row list
  indices = []  # csr column list
  data = []  # csr values
  objectives = []
  target = target.reset_index().set_index(['query_mode', 'query_size', 'query_index'])
  with open(topo_path, 'r') as f:
    for line in f:
      line = line.strip()
      splits = line.split(',')
      query_mode, query_size, query_index = splits[0].rsplit('/', 1)[1].rsplit("_", 3)[1:]
      query_config = (query_mode, int(query_size), int(query_index.split(".", 1)[0]))
      if query_config not in target.index:
        continue
      for kv in splits[1:]:
        k, v = kv.split(':')
        if k.startswith("clique") or k.startswith("cycle"):
          #if not k.startswith("clique"):
          continue
        idx = label_sequence_map.setdefault(k, len(label_sequence_map))
        indices.append(idx)
        data.append(int(v))
      indptr.append(len(indices))
      query_configs.append(query_config)
      objectives.append(target.loc[query_config][['enumerate_time', 'n_embeddings', 'plan_time']].values)
  sparse_features = csr_matrix((data, indices, indptr), dtype=int)
  objectives = np.array(objectives)
  #reg = Ridge(solver='sag', alpha=0.05, normalize=True).fit(sparse_features, objectives[:,1])
  print("dimensionality", len(label_sequence_map))
  for train_idx, validate_idx in KFold(n_splits=2, shuffle=True).split(sparse_features, objectives[:, 1]):
    reg = LinearRegression(normalize=True).fit(sparse_features[train_idx], objectives[:, 1][train_idx])
    print("Train R square =", reg.score(sparse_features[train_idx], objectives[:, 1][train_idx]))
    print("Validate R square =", reg.score(sparse_features[validate_idx], objectives[:, 1][validate_idx]))
    feature_weights = np.absolute(reg.coef_)
    top_coef_features = feature_weights.argsort()[:][::-1]
    feature_weight_rank = [None] * len(feature_weights)
    for i, feature_idx in enumerate(top_coef_features):
      feature_weight_rank[feature_idx] = i
    feature_id_name = [None] * len(label_sequence_map)
    for k, v in label_sequence_map.items():
      feature_id_name[v] = k
    print("----- top 20 feature weight -----")
    for idx in top_coef_features[:20]:
      print(feature_id_name[idx], reg.coef_[idx], sparse_features.getcol(idx).sum())
    feature_occurrence = np.array([sparse_features.getcol(i).sum() for i in range(0, len(label_sequence_map))])
    occurence_rank = feature_occurrence.argsort()[-1:-21:-1]
    print("----- top 20 feature occurence -----")
    for idx in occurence_rank:
      print(feature_id_name[idx], feature_weight_rank[idx], reg.coef_[idx], feature_occurrence[idx])


def topo_performance_analysis2(target, topo_path, baseline_file):
  topo = pd.read_csv(topo_path)
  topo_features = [x for x in topo.columns if x != "query"]
  query_config = topo['query'].str.rsplit("/", 1).str[1].str.rsplit("_", 3)
  topo["query_mode"] = query_config.str[1]
  topo["query_size"] = pd.to_numeric(query_config.str[2])
  topo["query_index"] = pd.to_numeric(query_config.str[-1].str.split(".").str[0])
  topo.drop(columns=['query'], inplace=True)
  print("dense >>>>>>>>>>>>>>>>>>>>>>")
  print(topo.loc[topo['query_mode'] == "dense"].describe())
  print("sparse >>>>>>>>>>>>>>>>>>>>>>")
  print(topo.loc[topo['query_mode'] == "sparse"].describe())

  topo = topo.loc[topo['query_mode'] == "dense"]
  topo = topo.set_index(['query_mode', 'query_size', 'query_index']).sort_index()
  target = target.reset_index().set_index(['query_mode', 'query_size', 'query_index'])
  data = target.join(topo, how='inner')

  reg = LinearRegression().fit(data[topo_features], data['enumerate_time'])
  print("R square =", reg.score(data[topo_features], data['enumerate_time']))
  print("model", reg.coef_, reg.intercept_)
  clf = MLPRegressor(solver='lbfgs', learning_rate='adaptive', random_state=None, max_iter=1000).fit(data[topo_features].values,
                                                                                                     data['n_embeddings'].values)
  print("R-square", clf.score(data[topo_features].values, data['n_embeddings'].values))

  if baseline_file is not None:
    baseline = read(baseline_file)[["enumerate_time"]].dropna().reset_index().set_index(['query_mode', 'query_size', 'query_index'])
    data = baseline.join(data, lsuffix='_target', rsuffix='_base', how='inner')

    data["time_diff"] = data['enumerate_time_target'] - data['enumerate_time_base']
    clf = MLPRegressor(solver='lbfgs', learning_rate='adaptive', random_state=None, max_iter=1000).fit(data[topo_features].values,
                                                                                                       data['time_diff'].values)
    print("time diff R-square", clf.score(data[topo_features].values, data['time_diff'].values))
    print("======= largest time diff =======")
    print(data.nlargest(10, 'time_diff')[['time_diff', *topo_features]])
    print("======= smallest time diff =======")
    print(data.nsmallest(10, 'time_diff')[['time_diff', *topo_features]])
    print("======= largest count =======")
    print(data.nlargest(10, 'n_embeddings')[['n_embeddings', *topo_features]])
    print("======= largest time =======")
    print(data.nlargest(10, 'enumerate_time_target')[['enumerate_time_target', *topo_features]])
    print("======= smallest time =======")
    print(data.nsmallest(10, 'enumerate_time_target')[['enumerate_time_target', *topo_features]])
    #reg = LinearRegression().fit(data[topo_features], data['enumerate_time_target'] - data['enumerate_time_base'])
    #print("R square =", reg.score(data[topo_features], data['enumerate_time_target'] - data['enumerate_time_base']))
    #print("model", reg.coef_, reg.intercept_)


def analyze_profile(folder, corr=True, candidate=False, cost_time_reg=True):
  profiles = [name for name in listdir(folder) if name != "log" and not name.endswith(".swp") and not name.endswith(".log")]

  df = None
  for name in tqdm(profiles):
    df = read_profile(osp.join(folder, name), df)
  #print(df)
  pd.set_option('display.max_rows', None)  # print all rows without ellipsis
  if corr:
    print("================================== Correlation ==================================")
    print(df[['cost', 'si_count', 'si_input', 'time']].corr())
    print("=============================== Correlation By Op ===============================")
    print(df.groupby('op')[['cost', 'si_count', 'si_input', 'time']].corr(min_periods=1))
    print("========================== Correlation By Op w/ Parent ==========================")
    print(df.groupby(['op', 'n_parents'])[['cost', 'si_count', 'si_input', 'time']].corr(min_periods=1))
    print("============================== Correlation By Step ==============================")
    print(df.groupby('step')[['cost', 'si_count', 'si_input', 'time']].corr(min_periods=1))
    print("=========================== Correlation By Query Mode ===========================")
    print(df.groupby('mode')[['cost', 'si_count', 'si_input', 'time']].corr(min_periods=1))
    print("========================== Correlation By Parent Count ==========================")
    print(df.groupby('n_parents')[['cost', 'si_count', 'si_input', 'time']].corr(min_periods=1))
  if candidate:
    print("========================== Candidate si effect ==========================")
    print(df['candidate_si_effect'].describe())
    print(df.groupby('step')['candidate_si_effect'].describe())
    print(df.groupby('op')['candidate_si_effect'].describe())
  if cost_time_reg:
    reg = LinearRegression().fit(df[['cost']], df['time'])
    print("R square =", reg.score(df[['cost']], df['time']))
    print("model", reg.coef_, reg.intercept_)
    print(df[['cost', 'time']].loc[df['time'] > 1000].nsmallest(10, 'cost'))
    print(df[['cost', 'time']].loc[df['time'] > 1000].nlargest(10, 'cost'))


def plan_time_analysis(target, baseline_file):
  compare_time(target, baseline_file, ['plan_time', 'enumerate_time'])
  print(target.nlargest(10, 'plan_time'))
  print("< 1 second query total plan time", target['plan_time'].loc[target['enumerate_time'] < 1].sum(), "total enumerate time", target['enumerate_time'].loc[target['enumerate_time'] < 1].sum(), "count",len(target['plan_time'].loc[target['enumerate_time'] < 1]))
  print("total plan time", target['plan_time'].sum(), "total enumerate time", target['enumerate_time'].sum(), "count",len(target['plan_time']))


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
    compare_time(target, args.time)
  if args.missing:
    print("Total elapsed backtracking time", target['elapsed_time'].sum())
    get_missing_query_config(target, args.missing)
  if args.config:
    get_config(target, args.config)
  if args.parallel:
    parallel_analysis(target)
  if args.topo:
    topo_performance_analysis(target, args.topo, args.time)
  if args.plan:
    plan_time_analysis(target, args.plan)

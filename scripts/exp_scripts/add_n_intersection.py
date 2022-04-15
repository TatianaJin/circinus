import argparse
import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import sys

from os import path as osp
from os import listdir
from sklearn.linear_model import LinearRegression, ElasticNet, Ridge
from sklearn.model_selection import KFold
from sklearn.neural_network import MLPRegressor
from scipy.sparse import csr_matrix
from sys import stderr
from tqdm import tqdm



def read(path, query_mode = "all"):
  f = pd.read_csv(path, dtype={'query_size':'Int64', 'query_index':'Int64'})
  if query_mode != "all":
    f = f.loc[f['query_mode'] == query_mode]
  return f.set_index(['dataset', 'query_size', 'query_mode', 'query_index'])

def read_profile_log(prefix, query_size, algo, cover, df):
    for index, row in df.iterrows():
        n_intersections = 0
        min_n_intersections = 0
        profile_path = osp.join(prefix, "{0}_{1}_{2}_{3}_{4}_{5}_{4}_cc".format(index[0], index[2], query_size, index[3], algo, cover))
        with open(profile_path, 'r') as profile:
            for line in profile:
                s = line.split(' ')
                if s[0] == 'step_costs' or len(s) == 1 or s[1] == 'Plan' or s[1] == 'query':
                    continue
                elif s[1] == 'Partitions':
                    break
                ss = s[-1].split(',')
                n_intersections += int(ss[6])
        df.at[index, 'n_intersections'] = n_intersections


def add_n_intersection(df_file, size, algo, mode):

    df = read(df_file).dropna()
    read_profile_log(df_file + "_profile", size, algo, mode, df)
    df['n_intersections'] = df['n_intersections'].astype('uint64')
    #df = df.reset_index()
    df.to_csv("{0}_count".format(df_file))


if __name__ == '__main__':
    df_file = sys.argv[1]
    size = sys.argv[2]
    algo = sys.argv[3]
    mode = sys.argv[4]
    add_n_intersection(df_file, size, algo, mode)
    # dataset = "youtube2007"
    # for dataset in ["youtube2007", "human", "orkut", "friendster"]:
    #     for algo in ["cfl","gql"]:
    #         for size in [8,12]:
    #             for mode in ["dynamic"]:
    #                 for path in ["true", "false"]:
    #                     add_n_intersection(dataset, algo, size, mode, path)


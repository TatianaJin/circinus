#! /usr/bin/env python3

import argparse
from os import path as osp
import os
from sys import stderr

def get_args():
    parser = argparse.ArgumentParser(description='Benchmark')
    parser.add_argument('-p', '--project_dir', help='Circinus project dir.')
    parser.add_argument('-t', '--time_out', default=300, help='Query timeout.')
    parser.add_argument('-m', '--match', default='Benchmark', help='The executable for running circinus.')

    # opts for workload
    parser.add_argument('-b', '--batch', default=1024, help='Batch size.')
    parser.add_argument('-f', '--filter', default='cfl', help='Candidate filter.')
    parser.add_argument('-c', '--config', default='cfl_4.csv', help='The query configs.')
    parser.add_argument('-l', '--limit', default=0, help='Match limit. 0 means unlimited.')
    parser.add_argument('-s', '--strategy', default='dynamic', choices=['dynamic','eager','static','all'], help='Vertex cover strategy.')
    parser.add_argument('--profile', choices=['simple', 'si'])
    
    args = parser.parse_args()
    if args.project_dir is None:
        print("please set project_dir")
        parser.print_usage()
        exit(0)
    args.bin_dir=osp.join(args.project_dir,"build","tests")
    args.log_dir=osp.join(args.project_dir,"evaluation","log")
    args.config_dir=osp.join(args.project_dir,"evaluation","configs")
    args.run_log=osp.join(args.project_dir,"evaluation","run_log")
    return args,parser

def run_batch(args):
    config_path = osp.join(args.config_dir, args.config)
    if not osp.isfile(config_path):
        stderr.write("File not exists {0}\n".format(config_path))
        return
    log = osp.join(args.log_dir, "{0}_{1}_{2}_{3}_{4}_{5}{6}".format(args.match,args.limit,args.batch,args.filter,args.strategy,args.config,"" if args.profile is None else "_{0}".format(args.profile)))
    if args.profile is not None:
        log_dir = log
        os.mkdirs(log_dir, exist_ok=True)
        log = osp.join(log_dir, "log")
    with open(log, 'a') as log_f:
        log_f.write("dataset,query_size,query_mode,query_index,load_time,filter_time,plan_time,enumerate_time,n_embeddings,order\n")
    common_flags = "-match_limit {0} -batch_size {1} -filter {2} -vertex_cover {3} -output_file {4}".format(args.limit, args.batch, args.filter, args.strategy, log)
    with open(config_path, 'r') as config_f:
        for line in config_f:
            if line.startswith("dataset"):
                continue
            flags = line.split(',')
            query_id="_".join(flags[0:4])
            # TODO(byli): profile_flag
            profile_flag="" if args.profile is None else "-profile_file {0}".format(osp.join(log_dir,query_id))
            cmd = "timeout {0} {1} -dataset {2} -query_size {3} -query_mode {4} -query_index {5} -match_order '{6}' {7} {8}".format(args.time_out, osp.join(args.bin_dir,args.match), flags[0], flags[1], flags[2], flags[3], flags[4].strip(), profile_flag, common_flags)
            os.system(cmd)


""" default workloads """
def run_cfl(args):
    pass

if __name__ == '__main__':
    args,_= get_args()
    if not osp.exists((args.log_dir)):
        os.mkdir(args.log_dir)
    print(args)
    run_batch(args)

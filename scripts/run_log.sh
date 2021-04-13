#! /bin/bash

# SET YOUR OWN
proj_dir=/data/share/users/tati/subgraph_matching/circinus

time_limit=300

# default workload
batch_size=1024
match_limit=0
vertex_cover=dynamic
filter=cfl
config=cfl_4.csv

run_batch () {
  if [[ $# -gt 2 ]]; then config=$1;       fi
  if [[ $# -gt 2 ]]; then filter=$2;       fi
  if [[ $# -gt 2 ]]; then vertex_cover=$3; fi
  if [[ $# -gt 2 ]]; then match_limit=$4;  fi
  if [[ $# -gt 2 ]]; then batch_size=$5;   fi
  
  python3 run_benchmark.py -b $batch_size -f $filter -l $match_limit -s $vertex_cover -c $config -l $match_limit
}

run_cfl_4() {
  batch_size=1024
  run_batch cfl_4.csv cfl dynamic 0
  run_batch cfl_4.csv cfl static 0
  run_batch cfl_4.csv cfl all 0
}

run_cfl_8() {
  batch_size=1024
  run_batch cfl_8.csv cfl dynamic 0
  run_batch cfl_8.csv cfl static 0
  run_batch cfl_8.csv cfl all 0
}

$@

#! /bin/bash

# SET YOUR OWN
proj_dir=/data/share/users/tati/subgraph_matching/circinus

time_limit=300

# default workload
batch_size=1024
match_limit=0
vertex_cover=dynamic
filter=cfl
order=cfl
config=cfl_4.csv
partition=1
upg=1
ipp=1
pqv=none

run_batch () {
  if [[ $# -gt 0 ]]; then config=$1;       fi
  if [[ $# -gt 1 ]]; then filter=$2;       fi
  if [[ $# -gt 2 ]]; then order=$3;        fi
  if [[ $# -gt 3 ]]; then vertex_cover=$4; fi
  if [[ $# -gt 4 ]]; then match_limit=$5;  fi
  if [[ $# -gt 5 ]]; then batch_size=$6;   fi
  
  python3 run_benchmark.py -p $proj_dir -t $time_limit -c $config \
    -b $batch_size -f $filter -o $order -s $vertex_cover -l $match_limit \
    --partition $partition --upg $upg --ipp $ipp --pqv $pqv
}

run_cfl_4() {
  batch_size=1024
  run_batch cfl_4.csv cfl cfl dynamic 0
  run_batch cfl_4.csv cfl cfl static 0
  run_batch cfl_4.csv cfl cfl all 0
}

run_cfl_8() {
  batch_size=1024
  run_batch cfl_8.csv cfl cfl dynamic 0
  run_batch cfl_8.csv cfl cfl static 0
  run_batch cfl_8.csv cfl cfl all 0
}

peu2005() {
  if [[ $# -gt 0 ]]; then partition=$1; else partition=20; fi
  run_batch eu2005_8_all.csv
}

$@

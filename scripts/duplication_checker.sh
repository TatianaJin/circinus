#!/bin/bash
# check all query graphs in $path using $checker, and log to $logdir
path="${1:-/data/share/users/qlma/query-graph-output}"
checker="${2:-/home/qlma/circinus/build/tests/Benchmark}"
logdir="${3:-/home/qlma/log/}"
let cnt=0
for dir1 in $path/*; do
  for dir2 in $path/*; do
    file1=$(basename $dir1)
    file2=$(basename $dir2)
    if [ $file1 = $file2 ]
    then
      break
    else
      prefix1=${file1%_*}
      prefix2=${file2%_*}
      if [ $prefix1 = $prefix2 ]
      then
        $checker -naive_run=1 -naive_datagraph=$dir1 -naive_querygraph=$dir2 -log_dir=$logdir | grep 'MATCH!!!'
        if [ $? -eq 0 ]
        then
          echo $dir1
          break
        fi
      fi
    fi
  done
  if [ $((cnt%10)) -eq 0 ]
  then
    echo $cnt
  fi
  let cnt=cnt+1
done


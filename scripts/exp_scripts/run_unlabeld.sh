#!/usr/bin/env bash

project_dir=$CIRCINUS_HOME # please set project dir
output_file=${project_dir}/exp_results_artifact/exp/data/compare_unlabeled/
data_dir=${project_dir}/subgraph_matching_datasets
vertex_cover=dynamic
num_cores=32
batch_size=128
bin=${project_dir}/build/tests/Benchmark

batch_run() {
	echo dataset,elapsed_execution_time,filter_time,plan_time,enumerate_time,n_embeddings,order,max_task_time >> ${output_file}/batch_${dataset}_${vertex_cover}_${num_cores}${log}
	for query_index in {1..8}; do
		$bin --data_dir ${data_dir} --dataset ${dataset} --query_index ${query_index} --match_limit 0 --num_cores ${num_cores} --output_file ${output_file}/batch_${dataset}_${num_cores}${log} --verbosity 0 --batch_size ${batch_size} --profile 1 --unlabeled 1 --break_symmetry 1
	done
}

dataset=${1}
log=_unlabeled
batch_run

project_dir=$CIRCINUS_HOME
output_file=${project_dir}/exp_results_artifact/exp/data/parallelization/
data_dir=${project_dir}/subgraph_matching_datasets
match_order=cfl
filter=cfl
pqv=cc
num_cores=${1}
batch_size=128
bin=${project_dir}/build/tests/Benchmark

batch_run() {
	echo dataset,query_size,query_mode,query_index,elapsed_execution_time,filter_time,plan_time,enumerate_time,n_embeddings,order,max_task_time >> ${output_file}/batch_${dataset}_${size}_${match_order}_${filter}_${partition}_${vertex_cover}_${pqv}_${num_cores}${log}

	echo ${bin} --batch_run 1 --data_dir ${data_dir} --dataset ${dataset} --query_size ${size} --match_order ${match_order} --filter ${filter} --partition ${partition} --match_limit 0 --vertex_cover ${vertex_cover} --pqv ${pqv} --batch_run 1 --num_cores ${num_cores} --upg 1 --output_file ${output_file}/batch_${dataset}_${size}_${match_order}_${filter}_${partition}_${vertex_cover}_${pqv}_${num_cores} --verbosity 0 --batch_size ${batch_size} --profile 1
	$bin --break_symmetry 1 --batch_run 1 --data_dir ${data_dir} --dataset ${dataset} --query_size ${size} --query_mode ${query_mode} --query_index ${query_index} --match_order ${match_order} --filter ${filter} --partition ${partition} --match_limit 0 --vertex_cover ${vertex_cover} --pqv ${pqv} --num_cores ${num_cores} --upg 1 --output_file ${output_file}/batch_${dataset}_${size}_${match_order}_${filter}_${partition}_${vertex_cover}_${pqv}_${num_cores}${log} --verbosity 0 --batch_size ${batch_size} --profile 1
}

dataset=friendster
partition=10
size=16
vertex_cover=dynamic
log=_parallel
batch_run

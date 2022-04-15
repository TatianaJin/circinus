project_dir=$CIRCINUS_HOME
output_file=${project_dir}/exp_results_artifact/exp/data/vc_influence/
data_dir=${project_dir}/subgraph_matching_datasets
match_order=online
filter=online
partition=20
pqv=cc
num_cores=32
batch_size=128
bin=${project_dir}/build/tests/Benchmark

batch_run() {
	echo $bin --dataset ${dataset} --query_size ${size} --match_order ${match_order} --filter ${filter} --partition ${partition} --match_limit 0 --vertex_cover ${vertex_cover} --pqv ${pqv} --batch_run 1 --num_cores ${num_cores} --upg 1 --output_file ${output_file}/batch_${dataset}_${size}_${match_order}_${filter}_${partition}_${vertex_cover}_${pqv}_${num_cores}${log} --verbosity 0 --batch_size ${batch_size} --profile 1
	$bin --data_dir ${data_dir} --dataset ${dataset} --query_size ${size} --match_order ${match_order} --filter ${filter} --partition ${partition} --match_limit 0 --vertex_cover ${vertex_cover} --pqv ${pqv} --batch_run 1 --num_cores ${num_cores} --upg 1 --output_file ${output_file}/batch_${dataset}_${size}_${match_order}_${filter}_${partition}_${vertex_cover}_${pqv}_${num_cores}${log} --verbosity 0 --batch_size ${batch_size} --profile 1 --path_card=${path_card} --break_symmetry 1 --time_limit 1801
}

batch_run_mode() {
	echo dataset,query_size,query_mode,query_index,elapsed_execution_time,filter_time,plan_time,enumerate_time,n_embeddings,order,max_task_time >> ${output_file}/batch_${dataset}_${size}_${match_order}_${filter}_${partition}_${vertex_cover}_${pqv}_${num_cores}${log}
	for query_mode in {dense,sparse}; do
		for i in {1..100}; do
			query_index=${i}
			timeout 1801 $bin --data_dir ${data_dir} --dataset ${dataset} --query_size ${size} --query_mode ${query_mode} --query_index ${query_index} --match_order ${match_order} --filter ${filter} --partition ${partition} --match_limit 0 --vertex_cover ${vertex_cover} --pqv ${pqv} --num_cores ${num_cores} --upg 1 --output_file ${output_file}/batch_${dataset}_${size}_${match_order}_${filter}_${partition}_${vertex_cover}_${pqv}_${num_cores}${log} --verbosity 0 --batch_size ${batch_size} --profile 1 --break_symmetry 1 --path_card=${path_card}
		done
	done
}

vertex_cover=${1}
size=${2}
dataset=${3}
log=_dynamic_sensitive

if [[ "$dataset" == "human" ]]
then
	partition=3
	path_card=false
	batch_run_mode
fi

if [[ "$dataset" == "youtube2007" ]]
then
	partition=5
	path_card=false
	batch_run_mode
fi


if [[ "$dataset" == "orkut" ]]
then
	partition=3
	path_card=false
	batch_run
fi


if [[ "$dataset" ==  "friendster" ]]
then
	partition=10
	path_card=false
	batch_run
fi


# 200 queries of each query vertex size
# Performance comparison
# labeled query
# Figure6
single_thread_compare(){
	# compare with cfl/gql alogrithm
	for size in {8,12,16}; do
			bash run_human_single.sh cfl ${size} ${i} # cfl is Circinus-CFL
			bash run_human_single.sh gql ${size} ${i} # gql is Circinus-GQL
			bash run_human_single.sh online ${size} ${i} # online is Circinus
	done

	for size in {8,12,16}; do
		bash run_youtube_single.sh cfl ${size} ${i} # cfl is Circinus-CFL
		bash run_youtube_single.sh gql ${size} ${i} # gql is Circinus-GQL
		bash run_youtube_single.sh online ${size} ${i} # online is Circinus
	done
}

mutli_thread_compare(){ # 32 hype-threads
	# compare with peregrine
	for size in {6,8}; do
		for algo in {cfl,online}; do
			bash run_youtube_multi.sh ${algo} ${size}
			bash run_orkut_multi.sh ${algo} ${size}
		done
	done
}

# unlabeled query
# Figure7


# Effectiveness of Redundancy Reduction
# Figure8 uses the result in single_thread_compare
# Figure9
redundancy_reduction(){
	# orkut
	for algo in {cfl,gql}; do
		for size in {8,12,16}; do
			for vc in {dynamic,none}; do # dynamic is Circinus, none is Circinus-
				bash run_orkut_redundancy_reduction.sh ${algo} ${size} ${vc}
			done
		done
	done

	# youtube
	for algo in {cfl,gql}; do
		for vc in {dynamic,none}; do # dynamic is Circinus, none is Circinus-
			bash run_youtube_redundancy_reduction.sh ${algo} 8 ${vc}
		done
	done

	# friendster
	for algo in {cfl,gql}; do
		for vc in {dynamic,none}; do # dynamic is Circinus, none is Circinus-
			bash run_friendster_redundancy_reduction.sh ${algo} 8 ${vc}
		done
	done
}

# Influence of Vertex Covers
# Figure10
vc_influence() {
	for dataset in {human,youtube2007,orkut,friendster}; do
		for vc in {dynamic,static}; do # dynamic is Circinus, static is MWVC
			for size in {8,12}; do
				bash run_vc_influence.sh ${vc} ${size} ${dataset}
			done
		done
	done
}

# cardinality_sensitivity
# Figure11
cardinality_sensitivity(){
	for dataset in {human,youtube2007,orkut,friendster}; do
		for path_card in {true,false}; do # true is Path, false is Product
			for size in {8,12}; do
				for algo in {cfl,gql}; do
					bash run_vc_influence.sh ${algo} ${path_card} ${dataset} ${size}
				done
			done
		done
	done
}

# Performance of Parallelization
# Figure12
parallelization() {
	for num_cores in {1,2,4,8,16,32}; do
		bash run_parallelization.sh ${num_cores}
	done
}

single_thread_compare
mutli_thread_compare
redundancy_reduction
vc_influence
cardinality_sensitivity
parallelization


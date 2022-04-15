#./process_log.py plot -m sum exp/data/compare_single_thread/human_youtube2007_online_single.csv 'exp/data/compare_single_thread/batch_human_*_cfl_cfl_3_dynamic_cc_1_single' 'exp/data/compare_single_thread/batch_youtube2007_*_cfl_cfl_5_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/batch_human_*_gql_gql_3_dynamic_cc_1_single' 'exp/data/compare_single_thread/batch_youtube2007_*_gql_gql_5_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/CFL.*.*.log.single' 'exp/data/compare_single_thread/GQLfs.*.*.log.single' -t Circinus Circinus-CFL Circinus-CFL Circinus-GQL Circinus-GQL CFL GQL -o exp/labeled_single_thread_sum_time
#./process_log.py plot -m sum exp/data/compare_peregrine/peregrine.orkut.6.log exp/data/compare_peregrine/peregrine.orkut.8.log exp/data/compare_peregrine/peregrine.youtube2007.6.log exp/data/compare_peregrine/peregrine.youtube2007.8.log exp/data/compare_peregrine/batch_orkut_6_cfl_cfl_3_dynamic_cc_32_multi exp/data/compare_peregrine/batch_orkut_8_cfl_cfl_3_dynamic_cc_32_multi exp/data/compare_peregrine/batch_youtube2007_6_cfl_cfl_5_dynamic_cc_32_multi exp/data/compare_peregrine/batch_youtube2007_8_cfl_cfl_5_dynamic_cc_32_multi exp/data/compare_peregrine/batch_orkut_6_online_online_3_dynamic_cc_32_multi exp/data/compare_peregrine/batch_orkut_8_online_online_3_dynamic_cc_32_multi exp/data/compare_peregrine/batch_youtube2007_6_online_online_5_dynamic_cc_32_multi exp/data/compare_peregrine/batch_youtube2007_8_online_online_5_dynamic_cc_32_multi -t Peregrine Peregrine Peregrine Peregrine Circinus-CFL Circinus-CFL Circinus-CFL Circinus-CFL Circinus Circinus Circinus Circinus -o exp/labeled_multi_thread_sum_time
#./process_log.py plot -m si exp/data/compare_single_thread/human_youtube2007_online_single.csv 'exp/data/compare_single_thread/batch_human_*_cfl_cfl_3_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/batch_youtube2007_*_cfl_cfl_5_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/batch_human_*_gql_gql_3_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/batch_youtube2007_*_gql_gql_5_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/CFL.*.*.log.single' 'exp/data/compare_single_thread/GQLfs.*.*.log.single' -t Circinus Circinus-CFL Circinus-CFL Circinus-GQL Circinus-GQL CFL GQL -o exp/redundancy_verification_intersection
#./process_log.py plot -m time exp/data/compare_single_thread/human_youtube2007_online_single.csv 'exp/data/compare_single_thread/batch_human_*_cfl_cfl_3_dynamic_cc_1_single' 'exp/data/compare_single_thread/batch_youtube2007_*_cfl_cfl_5_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/batch_human_*_gql_gql_3_dynamic_cc_1_single' 'exp/data/compare_single_thread/batch_youtube2007_*_gql_gql_5_dynamic_cc_1_single_count' 'exp/data/compare_single_thread/CFL.*.*.log.single' 'exp/data/compare_single_thread/GQLfs.*.*.log.single' -t Circinus Circinus-CFL Circinus-CFL Circinus-GQL Circinus-GQL CFL GQL -o exp/redundancy_verification_time
#./process_log.py plot -m elapsed 'exp/data/cardinality_sensitivity/*cfl*true*' 'exp/data/cardinality_sensitivity/*gql*true*' 'exp/data/cardinality_sensitivity/*cfl*false*' 'exp/data/cardinality_sensitivity/*gql*false*' -t CFL-Path GQL-Path CFL-Product GQL-Product -o exp/card_method_exp_all_backtracking_time --legend_pos='lower right'
#./process_log.py plot -m si 'exp/data/cardinality_sensitivity/*cfl*true*' 'exp/data/cardinality_sensitivity/*gql*true*' 'exp/data/cardinality_sensitivity/*cfl*false*' 'exp/data/cardinality_sensitivity/*gql*false*' -t CFL-Path GQL-Path CFL-Product GQL-Product -o exp/card_method_exp_all_intersection
#./process_log.py plot -m elapsed 'exp/data/vc_incluence/batch_*dynamic_cc*' 'exp/data/vc_incluence/batch_*static_cc*' -t Circinus MWVC --ylim 5e-4 1e3 -o exp/vc_incluence_backtracking_time --figsize 14 6
#./process_log.py plot_uncompressed -m elapsed 'exp/data/redundancy_multithread/*8*' --ylim 5e-4 1e3 -o exp/redundancy_8_time
#./process_log.py plot_uncompressed -m si 'exp/data/redundancy_multithread/*8*' --ylim 5e3 1e10 -o exp/redundancy_8_count
#./process_log.py plot_uncompressed -m elapsed 'exp/data/redundancy_multithread/*orkut*' --ylim 5e-4 1e3 -o exp/redundancy_orkut_time -x query
#./process_log.py plot_uncompressed -m si 'exp/data/redundancy_multithread/*orkut*' --ylim 5e3 1e10 -o exp/redundancy_orkut_count -x query
./process_log.py motivation exp/data/motivation/unlabel_motivation exp/data/motivation/motivation_labeled.csv -o motivation

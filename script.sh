
#  DIR STATS COUNTS PRIOR
cd analiza
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
do
	python -m generator --bins 2000 --clusters 200 --tree_size 19 --snv_in_cn_proportion 0.5 --cluster_size 1 --snv_events_per_edge 40 --sequencing_error 0.00001 --coverage 10 --read_success_prob 0.63

  docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.01 --real_breakpoints 1 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_real_brkp_1000_0.01"



  docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.1 --real_breakpoints 1 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_real_brkp_1000_0.1"


	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 10000 --counts_penalty_s1 10000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.01 --real_breakpoints 1 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_real_brkp_10000_0.01"



  docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 10000 --counts_penalty_s1 10000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.1 --real_breakpoints 1 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_real_brkp_10000_0.1"







  docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.01 --real_breakpoints 0 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_unknown_brkp_1000_0.01"



  docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.1 --real_breakpoints 0 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_unknown_brkp_1000_0.1"


	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 10000 --counts_penalty_s1 10000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.01 --real_breakpoints 0 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_unknown_brkp_10000_0.01"



  docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester clustered.py  --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 10000 --counts_penalty_s1 10000 --sequencing_error 0.00001 --pt_inf_iters 1000 --param_inf_iters 500000 --tree_structure_prior_k1 0.1 --real_breakpoints 0 --snv_clustered 0
	docker run -v /home/neuro/Desktop/snv-calling/analiza:/data tc360950/conset/tester tester.py --cbs_min_cells 3  --snv_constant 0.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.00001 --pt_inf_iters 1000000 --param_inf_iters 500000 --tree_structure_prior_k1 1 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../super_results2 --prefix "unknown_complex_inf_attachment_full_genotypes_unknown_brkp_10000_0.1"

done

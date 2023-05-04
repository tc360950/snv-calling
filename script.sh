
#  DIR STATS COUNTS PRIOR
cd $1
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
do
	python -m generator --bins 2000 --clusters 200 --tree_size 19 --snv_in_cn_proportion 0.5 --cluster_size 1 --snv_events_per_edge 40 --sequencing_error 0.00001 --coverage 10 --read_success_prob 0.63

	docker run -v /home/tc360950/new_sims/snv-calling/$1:/data tc360950/conset:1.0.0 clustered.py --cbs_min_cells 3  --snv_constant 1.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 $3 --counts_penalty_s1 $3 --sequencing_error 0.0001 --pt_inf_iters 100000 --param_inf_iters 100000 --tree_structure_prior_k1 $4 --real_breakpoints 1 --snv_clustered 0
	python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix "not_clustered_${3}_${4}"


  docker run -v /home/tc360950/new_sims/snv-calling/$1:/data tc360950/conset:1.0.0 clustered.py --cbs_min_cells 3  --snv_constant 1.0 --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 $3 --counts_penalty_s1 $3 --sequencing_error 0.0001 --pt_inf_iters 100000 --param_inf_iters 100000 --tree_structure_prior_k1 $4 --real_breakpoints 1 --snv_clustered 1
	python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix "not_clustered_${3}_${4}"

done

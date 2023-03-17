
cd $1
python -m generator --bins 10000 --clusters 1300 --tree_size 19 --snv_in_cn_proportion 0.5 --cluster_size 1 --snv_events_per_edge 40 --sequencing_error 0.00001 --coverage 0.077 --read_success_prob 0.63

docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.8 unclustered.py --seed 123  --cbs_min_cells 3 --min_coverage 5  --max_coverage 15 --estimate_snv_constant True --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 5000 --param_inf_iters 5000 --recalculate_cbs True --tree_structure_prior_k1 1.0

python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix unclustered_5_True

docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.8 unclustered.py --seed 123  --cbs_min_cells 3 --min_coverage 5  --max_coverage 15 --estimate_snv_constant True --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 5000 --param_inf_iters 5000  --tree_structure_prior_k1 1.0

python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix unclustered_5_True

docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.8 unclustered.py --seed 123  --cbs_min_cells 3 --min_coverage 5  --max_coverage 15 --estimate_snv_constant True --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 5000 --param_inf_iters 5000 --tree_structure_prior_k1 1.0

python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix unclustered_5_True


docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.8 unclustered.py --seed 4  --cbs_min_cells 3 --min_coverage 5  --max_coverage 15 --estimate_snv_constant True --snv_candidates 100 --output_dir /data/out --data_dir /data  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 5000 --param_inf_iters 5000  --tree_structure_prior_k1 1.0

python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix unclustered_5_True

#	docker run -v /home/tc360950/new_sims/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.8 unclustered.py --seed $i  --cbs_min_cells 3 --min_coverage 3 --max_coverage 12 --snv_candidates 100 --output_dir /data/out --data_dir /data --estimate_snv_constant True --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 100000 --param_inf_iters 100000 --tree_structure_prior_k1 1.0
#
#	python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix unclustered_3_True
#
#	docker run -v /home/tc360950/new_sims/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.8 unclustered.py --seed $i  --cbs_min_cells 3 --min_coverage 3 --max_coverage 12 --snv_candidates 100 --output_dir /data/out --data_dir /data --snv_constant 0.0  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 100000 --param_inf_iters 100000 --tree_structure_prior_k1 1.0
#
#	python -m generator --simulation_dir ./ --stats_dir ../$2 --prefix unclustered_3_False

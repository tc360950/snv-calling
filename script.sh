#export DOCKER_HOST=unix://$XDG_RUNTIME_DIR/docker.sock


#  DIR STATS
for i in 1
do
cd $1
	python -m generator --bins 10000 --clusters 1000 --tree_size 20 --snv_in_cn_proportion 0.5 --random_attachment False --cluster_size 1 --snv_events_per_edge 40 --sequencing_error 0.00001 --coverage 0.077 --read_success_prob 0.63

	docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.6 unclustered.py --seed $i  --cbs_min_cells 3 --min_coverage 5 --snv_candidates 100 --output_dir /data/out --data_dir /data --estimate_snv_constant True  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 10000 --param_inf_iters 10000

	python -m generator --simulation_dir ./ --stats_dir ../$2

	docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.6 unclustered.py --seed $i  --cbs_min_cells 3 --min_coverage 5 --snv_candidates 100 --output_dir /data/out --data_dir /data --snv_constant 0.0  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 10000 --param_inf_iters 10000

	python -m generator --simulation_dir ./ --stats_dir ../$2


	docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.6 unclustered.py --seed $i  --cbs_min_cells 3 --min_coverage 10 --snv_candidates 100 --output_dir /data/out --data_dir /data --estimate_snv_constant True  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 10000 --param_inf_iters 10000

	python -m generator --simulation_dir ./ --stats_dir ../$2

	docker run -v /home/neuro/Desktop/snv-calling/$1:/data tc360950/conet_snv:p_inference_v1.6 unclustered.py --seed $i  --cbs_min_cells 3 --min_coverage 10 --snv_candidates 100 --output_dir /data/out --data_dir /data --snv_constant 0.0  --counts_penalty_s2 1000 --counts_penalty_s1 1000 --sequencing_error 0.0001 --pt_inf_iters 10000 --param_inf_iters 10000

	python -m generator --simulation_dir ./ --stats_dir ../$2
done

 export DOCKER_HOST=unix://$XDG_RUNTIME_DIR/docker.sock


# SNVCONSTANT DIR STATS
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
do
cd $2
	python -m generator --bins 2000 --clusters 80 --tree_size 40 --snv_in_cn_proportion 0.5 --random_attachment False --cluster_size 100 --snv_events_per_edge 25


	docker run -v /home/tc360950/new_sims/snv-calling/$2:/data tc360950/conet_snv_recov:latest --data_dir /data/ --output_dir /data/out --corrected_counts_file /data/cc --pt_inf_iters 100000 --param_inf_iters 50000 --counts_penalty_s1 10000 --counts_penalty_s2 10000  --seed $i --tries 10

	python -m generator --simulation_dir ./ --stats_dir ../$3 --postfix 0
	python -m generator --simulation_dir ./ --stats_dir ../$3 --postfix 1
	python -m generator --simulation_dir ./ --stats_dir ../$3 --postfix 2
	python -m generator --simulation_dir ./ --stats_dir ../$3 --postfix 3
	python -m generator --simulation_dir ./ --stats_dir ../$3 --postfix 4
	python -m generator --simulation_dir ./ --stats_dir ../$3 --postfix 5
done

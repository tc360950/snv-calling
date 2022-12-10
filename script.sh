#export DOCKER_HOST=unix://$XDG_RUNTIME_DIR/docker.sock


# SNVCONSTANT DIR STATS
for i in 1
do
cd $2
	python -m generator --bins 50 --clusters 40 --tree_size 20 --snv_in_cn_proportion 0.5 --random_attachment False --cluster_size 100 --snv_events_per_edge 2

	docker run -v /home/neuro/Desktop/snv-calling/$2:/data tc360950/conet_snv:v2 --data_dir /data/ --output_dir /data/out --corrected_counts_file /data/cc --pt_inf_iters 100000 --param_inf_iters 30000 --counts_penalty_s1 10000 --counts_penalty_s2 10000  --snv_constant $1 --seed $i --tries 2

	python -m generator --simulation_dir ./ --stats_dir ../$3
done

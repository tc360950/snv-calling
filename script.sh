systemctl --user enable docker
export DOCKER_HOST=unix://$XDG_RUNTIME_DIR/docker.sock



for i in 1 2 3 4 5
do
	echo "Start $i"

	python -m generator --bins 100 --clusters 20 --tree_size 10 --snv_in_cn_proportion 0.5 --random_attachment False --cluster_size 100

	docker run -v /home/tc360950/new_sims/snv-calling/${SIM_DIR}:/data tc360950/conet_snv:latest --data_dir /data/ \
		--output_dir /data/out --corrected_counts_file /data/cc --pt_inf_iters 50000 \
		--param_inf_iters 30000 --counts_penalty_s1 10000 --counts_penalty_s2 10000 --snv_constant 0.0001 --seed $i

	python -m generator --simulation_dir ${SIM_DIR} --stats_dir ${STATS_DIR}
	echo "End $i"
done






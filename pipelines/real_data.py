import argparse
import pickle

import networkx as nx
import numpy as np
import pandas as pd
import solver.solver
from conet_integration.loader import CONETLoader
from core.model_data import CellsData
from prefect import Flow, task
from prefect.tasks.docker.containers import (
    CreateContainer,
    GetContainerLogs,
    StartContainer,
    WaitOnContainer,
)

DOCKER_HOST = ""

start = StartContainer(docker_server_url=DOCKER_HOST)
wait = WaitOnContainer(docker_server_url=DOCKER_HOST)
logs = GetContainerLogs()

parser = argparse.ArgumentParser(
    description="Run join CONET and SNV inference on real data."
)
parser.add_argument("--working_dir", type=str, default="/home/tomasz/Desktop/SNV/dump/")


@task
def attach_snvs_to_conet(dir: str) -> None:
    snv_to_bin = list(pd.read_csv(f"{dir}snv_to_bin", header=None)[0])
    loader = CONETLoader(dir=dir, snv_to_bin=snv_to_bin)
    tree_counts = loader.get_event_tree()
    attachment = loader.attachment
    d = np.transpose(pd.read_csv(f"{dir}d", header=0, sep=";").to_numpy())
    b = np.transpose(pd.read_csv(f"{dir}b", header=0, sep=";").to_numpy())

    cells = CellsData(
        d=d.astype(int),
        b=b.astype(int),
        attachment=[attachment[i] for i in range(0, len(attachment))],
        cell_cluster_sizes=[1 for _ in range(d.shape[0])],
    ).aggregate()

    inf_snv_tree = solver.solver.solve(cells, tree_counts)

    node_to_conv = {n1: n2 for n2, n1 in loader.node_to_conv.items()}
    inf_snv_tree.cn_event_tree = nx.relabel_nodes(
        inf_snv_tree.cn_event_tree, node_to_conv
    )
    inf_snv_tree.node_to_snvs = {
        node_to_conv[n]: snvs for n, snvs in inf_snv_tree.node_to_snvs.items()
    }
    with open("result", "wb") as f:
        pickle.dump(inf_snv_tree, f)


if __name__ == "__main__":
    args = parser.parse_args()

    with Flow("CONET SNV joint inference") as flow:
        create_conet_container = CreateContainer(
            image_name="tc360950/conet:latest",
            command=[
                "--data_dir=/data/",
                "--output_dir=/data/",
                "--corrected_counts_file=/data/cc",
            ],
            docker_server_url=DOCKER_HOST,
            volumes=["/data"],
            host_config={
                "binds": {
                    f"{args.working_dir}": {
                        "bind": "/data",
                        "mode": "rw",
                    }
                }
            },
        )

        conet_container = create_conet_container()
        start_conet = start(container_id=conet_container)
        wait_for_conet = wait(container_id=conet_container)
        conet_logs = logs(container_id=conet_container).set_upstream(wait_for_conet)

        snvs = attach_snvs_to_conet(args.working_dir).set_upstream(conet_logs)

    flow.run()

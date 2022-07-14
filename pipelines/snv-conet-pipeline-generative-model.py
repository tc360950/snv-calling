import argparse
import pickle

import numpy as np
from prefect import task, Flow, Parameter
from prefect.tasks.docker.containers import (
    CreateContainer,
    StartContainer,
    GetContainerLogs,
    WaitOnContainer,
)

from conet_integration.loader import CONETLoader
from conet_integration.utils import raw_cc_matrix_to_CONET_format
from core.model_data import CellsData
from generator.context import SNVGeneratorContext
from generator.gen import sample_cells_data, SNVModelGenerator
import solver.solver
from stats.statistics_calculator import StatisticsCalculator

DOCKER_HOST = ""

start = StartContainer(docker_server_url=DOCKER_HOST)
wait = WaitOnContainer(docker_server_url=DOCKER_HOST)
logs = GetContainerLogs()


@task
def generate_model_with_snv(working_dir: str, tree_size: int, clusters: int, bins: int, cluster_size: int) -> None:
    class GeneratorContext(SNVGeneratorContext):
        def number_of_bins(self):
            return bins

    model_generator = SNVModelGenerator(GeneratorContext())
    model = model_generator.generate_model(tree_size)

    while not any([np.min(n) == 0.0 for n in model.node_to_cn_profile.values()]):
        model = model_generator.generate_model(tree_size)

    cells_data, cc_matrix = sample_cells_data(clusters=clusters, cluster_size=cluster_size, model=model, ctxt=GeneratorContext())


    raw_cc_matrix_to_CONET_format(cc_matrix, model.tree).to_csv(f"{working_dir}cc", index=False)
    with open(f"{working_dir}cells_data", "w") as f:
        pickle.dump(cells_data, f)
    with open(f"{working_dir}event_tree", "wb") as f:
        pickle.dump(model.tree, f)


@task
def attach_snvs_to_conet(dir: str) -> None:
    with open(f"{dir}event_tree", "rb") as f:
        real_event_tree = pickle.load(f)
    with open(f"{dir}cells_data", "rb") as f:
        cells = pickle.load(f)
    tree_counts = CONETLoader(dir=dir).get_event_tree()

    attachment = CONETLoader(dir=dir).attachment

    cells = CellsData(d=cells.d.astype(int),b=cells.b.astype(int),
                      attachment=[attachment[i] for i in range(0, len(attachment))],
                      cell_cluster_sizes=cells.cell_cluster_sizes).aggregate()

    inf_snv_tree = solver.solver.solve(cells, tree_counts)

    print(StatisticsCalculator.calculate_stats(real_tree=real_event_tree, inferred_tree=inf_snv_tree))


with Flow("CONET SNV joint inference") as flow:
    cluster_size = Parameter('cluster_size')
    tree_size = Parameter('tree_size')
    clusters = Parameter('clusters')
    num_bins = Parameter('num_bins')
    working_dir = Parameter('working_dir')

    create_conet_container = \
        CreateContainer(image_name="tc360950/conet:latest",
                        command=["--data_dir=/data/", "--output_dir=/data/", "--corrected_counts_file=/data/cc"],
                        docker_server_url=DOCKER_HOST,
                        volumes=["/data"],
                        host_config={"binds": {
                            working_dir: {
                                'bind': '/data',
                                'mode': 'rw',
                            }
                        }})

    generate_snv_model = generate_model_with_snv(working_dir, tree_size, clusters, num_bins, cluster_size)
    conet_container = create_conet_container().set_upstream(generate_snv_model)
    start_conet = start(container_id=conet_container)
    wait_for_conet = wait(container_id=conet_container)
    conet_logs = logs(container_id=conet_container).set_upstream(wait_for_conet)
    attach_snvs_to_conet(working_dir).set_upstream(conet_logs)


parser = argparse.ArgumentParser(description='Run join CONET and SNV inference on generated data.')
parser.add_argument('--cluster_size', type=int)
parser.add_argument('--clusters', type=int)
parser.add_argument('--tree_size', type=int)
parser.add_argument('--num_bins', type=int)
parser.add_argument('--working_dir', type=str)

if __name__ == "__main__":
    args = parser.parse_args()

    flow.run(parameters=dict(cluster_size=args.cluster_size, clusters=args.clusters, tree_size=args.tree_size,num_bins=args.num_bins,working_dir=args.working_dir))

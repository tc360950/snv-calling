import json
from concurrent.futures import ProcessPoolExecutor
from typing import Tuple

import solver.solver
from generator.context import SNVGeneratorContext
from generator.model import SNVModel
from solver import CellsData

import numpy as np

from solver.input import EventTreeWithCounts
from stats.statistics_calculator import StatisticsCalculator

np.seterr(all='raise')

class GeneratorContext(SNVGeneratorContext):
   pass


def create_cn_matrix(cells: CellsData, model: SNVModel) -> np.ndarray:
    matrix = np.zeros(cells.d.shape)
    for cell in range(0, cells.d.shape[0]):
        matrix[cell,] = model.node_to_cn_profile[cells.attachment[cell]]
    return matrix.astype(int)


def simulate(a: Tuple[int ,int, int]):
    clusters = a[0]
    cluster_size = a[1]
    tree = a[2]
    CLUSTER_SIZE = cluster_size
    CLUSTERS = clusters
    iters = 5
    def generate_model():
        model = SNVModel(GeneratorContext())
        model.generate_structure(tree)
        cluster_to_cell_data = {}

        for c in range(0, CLUSTERS):
            print(f"Generating cluster num {c}")
            node = next(model.generate_cell_attachment())
            cell = model.generate_cell_in_node(node)
            for _ in range(0, CLUSTER_SIZE - 1):
                cell2 = model.generate_cell_in_node(node)
                cell.d = cell.d + cell2.d
                cell.b = cell.b + cell2.b
            cluster_to_cell_data[c] = cell

        cluster_to_size = dict([(c, CLUSTER_SIZE) for c in range(0, CLUSTERS)])
        attachment = dict([(c, cluster_to_cell_data[c].attachment) for c in range(0, CLUSTERS)])

        b = np.vstack([cluster_to_cell_data[c].b for c in range(0, CLUSTERS)])
        d = np.vstack([cluster_to_cell_data[c].d for c in range(0, CLUSTERS)])

        cells_data = CellsData(d=d, b=b, attachment=attachment, cell_cluster_sizes=cluster_to_size)
        return cells_data, model

    def generate_and_solve(*args, **kwargs):

        cells_data, model = generate_model()
        cn_profiles = create_cn_matrix(cells_data, model)
        while np.sum(cn_profiles == 0) == 0:
            cells_data, model = generate_model()
            cn_profiles = create_cn_matrix(cells_data, model)

        inferred_tree = solver.solver.solve(cells_data=cells_data,
                                            tree=EventTreeWithCounts(tree=model.event_tree.create_copy_without_snvs(), node_to_cn_profile=model.node_to_cn_profile))

        return StatisticsCalculator.calculate_stats(real_tree=model.event_tree, inferred_tree=inferred_tree)

    stats = [generate_and_solve() for _ in range(0, iters)]
    result = {}
    for key in stats[0].keys():
        result[key] = sum([s[key] for s in stats]) / iters
    return result, a

clusters = [50, 200, 500]
sizes = [10, 100]
tree = [10, 40]

a = []

for c in clusters:
    for s in sizes:
        for t in tree:
            a.append((c, s, t))

with open("results", "w") as f:
    with ProcessPoolExecutor(max_workers=10)as pool:
        for r in pool.map(simulate, a):
            print("OK")
            f.write(f"{r[1][0]};{r[1][1]};{r[1][2]};{json.dumps(r[0])}\n")
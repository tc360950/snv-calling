import json

from generator.context import SNVGeneratorContext
from generator.model import SNVModel
from solver import CellsData, MLSolver

import numpy as np


from stats.statistics_calculator import StatisticsCalculator


class GeneratorContext(SNVGeneratorContext):
   pass


def simulate(iters: int, clusters: int, cluster_size: int, tree: int):
    CLUSTER_SIZE = cluster_size
    CLUSTERS = clusters

    def generate():
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

        tree_with_snvs = model.event_tree
        model.event_tree = model.event_tree.create_copy_without_snvs()

        solver = MLSolver(cells_data, model, GeneratorContext())

        solver.insert_snv_events()

        return StatisticsCalculator.calculate_stats(model.event_tree, tree_with_snvs)

    stats = [generate() for _ in range(0, iters)]
    result = {}
    for key in stats[0].keys():
        result[key] = sum([s[key] for s in stats]) / iters
    return result

clusters = [50, 200, 500]
sizes = [10, 100]
tree = [10, 40]


with open("results", "w") as f :
    for c in clusters:
        for s in sizes:
            for t in tree:
                r = simulate(5, c, s, t)
                f.write(f"{c};{s};{t};{json.dumps(r)}\n")
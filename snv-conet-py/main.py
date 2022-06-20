import json

from generator.context import SNVGeneratorContext
from generator.model import SNVModel

import numpy as np


from stats.statistics_calculator import StatisticsCalculator


class GeneratorContext(SNVGeneratorContext):
   pass


model = SNVModel(GeneratorContext())
model.generate_structure(10)

CLUSTER_SIZE = 100
CLUSTERS = 10

from solver import CellsData, MLSolver

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

stats = StatisticsCalculator.calculate_stats(model.event_tree, tree_with_snvs)

print(json.dumps(stats))

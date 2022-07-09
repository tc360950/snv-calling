from core.event_tree import EventTree
from core.model_data import Parameters
from solver import CellsData, MLSolver
from solver.input import EventTreeWithCounts
from solver.parameters_inference import MMEstimator, NewtonRhapsonEstimator


def solve(cells_data: CellsData, tree: EventTreeWithCounts) -> EventTree:
    p: Parameters = MMEstimator.estimate(cells_data, tree)
    p = NewtonRhapsonEstimator(cells= cells_data, tree=tree).solve(p)
    return MLSolver(cd=cells_data, tree=tree, p=p).insert_snv_events([i for i in range(0, cells_data.d.shape[1])])

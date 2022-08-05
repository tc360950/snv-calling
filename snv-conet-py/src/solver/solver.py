from core.event_tree import EventTree
from core.logger import get_logger
from core.model_data import CellsData, EventTreeWithCounts, Parameters
from solver.ml_snv_solver import MLSolver
from solver.parameters_inference import MMEstimator, NewtonRhapsonEstimator

logger = get_logger(__name__)


def solve(cells_data: CellsData, tree: EventTreeWithCounts) -> EventTree:
    logger.info("Starting MM estimation of parameters...")
    p: Parameters = MMEstimator.estimate(cells_data, tree)
    logger.info(f"Estimated parameters: {p}. Starting Newton-Rhapson optimization...")
    p = NewtonRhapsonEstimator(cells=cells_data, tree=tree).solve(p)
    logger.info(f"Estimated parameters: {p}. Inferring SNV events...")
    return MLSolver(cd=cells_data, tree=tree, p=p).insert_snv_events(
        [i for i in range(0, cells_data.d.shape[1])]
    )

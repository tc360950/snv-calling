from typing import Callable, List, Optional, Set, Tuple

from core.event_tree import EventTree
from core.logger import get_logger
from core.model_data import CellsData, EventTreeWithCounts, Parameters
from core.types import CNEvent, EventTreeRoot, SNVEvent
from solver.likelihood_calculator import LikelihoodCalculator

logger = get_logger(__name__)


class MLSolver:
    cell_data: CellsData
    tree: EventTreeWithCounts
    calculator: LikelihoodCalculator

    def __init__(self, cd: CellsData, tree: EventTreeWithCounts, p: Parameters):
        self.cell_data = cd
        self.tree = tree
        self.calculator = LikelihoodCalculator(cell_data=cd, tree=tree, p=p)

    def insert_snv_events(self, snv_event_candidates: List[SNVEvent]) -> EventTree:
        calculator_refresh = self.calculator.get_likelihood_refresher()
        log_likelihood = calculator_refresh(None)
        logger.info(
            f"Starting SNV events inference, log-likelihood of tree with 0 SNVs: {log_likelihood}."
        )
        for snv in snv_event_candidates:
            snv_added = True
            while snv_added:  # One SNV can be added in multiple locations
                snv_added = False
                log_likelihood_snv = calculator_refresh(snv)
                best_node_for_snv, log_lik_with_snv = self.__choose_best_svn_location(
                    snv, calculator_refresh
                )
                if (
                    best_node_for_snv is not None
                    and log_lik_with_snv > log_likelihood_snv
                ):
                    logger.info(f"Adding SNV {snv} to the tree.")
                    self.tree.tree.add_snv_in_node(best_node_for_snv, snv)
                    calculator_refresh(snv)
                    snv_added = True

        logger.info(
            f"Finished SNV events inference. Log-likelihood {calculator_refresh(None)}."
        )
        return self.tree.tree

    def __choose_best_svn_location(
        self, snv: SNVEvent, calculator_refresh: Callable[[SNVEvent], float]
    ) -> Tuple[Optional[CNEvent], float]:
        logger.info(
            f"Choosing best location for SNV {snv}. "
            f"Likelihood without SNV: {calculator_refresh(snv)}."
        )
        location_to_lik = {}
        for node in self.__get_locations_for_snv(snv):
            self.tree.tree.add_snv_in_node(node, snv)
            location_to_lik[node] = calculator_refresh(snv)
            logger.info(f"Likelihood for node {node} is {location_to_lik[node]}...")
            self.tree.tree.delete_snv_from_node(node, snv)
            calculator_refresh(snv)

        best_node, best_lik = next(
            iter(
                [
                    (n, lik)
                    for n, lik in location_to_lik.items()
                    if lik == max(location_to_lik.values())
                ]
            ),
            (None, 0.0),
        )
        logger.info(
            f"Best location is node {best_node} with log-likelihood {best_lik}."
        )
        return best_node, best_lik

    def __get_locations_for_snv(self, snv: SNVEvent) -> Set[CNEvent]:
        return {
            n
            for n in self.tree.tree.cn_event_tree.nodes
            if self.tree.tree.snv_can_be_added_in_node(snv, n) and n != EventTreeRoot
        }

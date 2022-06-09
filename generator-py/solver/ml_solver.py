from typing import Tuple

import generator
from generator import CNEvent, SNVEvent, EventTreeRoot, get_logger
from solver.cell_data import CellsData
from solver.likelihood_calculator import LikelihoodCalculator

logger = get_logger(__name__)


class MLSolver:
    cell_data: CellsData
    model: generator.SNVModel
    context: generator.SNVGeneratorContext
    calculator: LikelihoodCalculator

    def __init__(self, cd: CellsData, model: generator.SNVModel, context: generator.SNVGeneratorContext):
        self.cell_data = cd
        self.model = model
        self.context = context
        self.calculator = LikelihoodCalculator(cell_data=cd, model=model, context=context)

    def insert_snv_events(self) -> generator.SNVModel:
        log_likelihood = self.calculator.get_likelihood()
        logger.info(f"Starting SNV events inference, log-likelihood {log_likelihood}.")
        for snv in self.context.get_snv_event_candidates():
            best_node_for_snv, log_lik_with_snv = self.__choose_best_svn_location(snv)
            if log_lik_with_snv > log_likelihood:
                logger.info(f"Adding SNV {snv} to the tree.")
                self.model.event_tree.node_to_snvs.setdefault(best_node_for_snv, set()).add(snv)
                log_likelihood = log_lik_with_snv
        logger.info(f"Finished SNV events inference. Log-likelihood {log_likelihood}.")
        return self.model

    def __choose_best_svn_location(self, snv: SNVEvent) -> Tuple[CNEvent, float]:
        logger.info(f"Choosing best location for SNV {snv}...")
        location_to_lik = {}
        for node in self.model.event_tree.cn_event_tree.nodes:
            if node == EventTreeRoot:
                continue

            self.model.event_tree.node_to_snvs.setdefault(node, set()).add(snv)
            location_to_lik[node] = self.calculator.get_likelihood()
            self.model.event_tree.node_to_snvs[node].remove(snv)

        best_node = [n for n, l in location_to_lik.items() if l == max(location_to_lik.values())][0]
        logger.info(f"Best location is node {best_node} with log-likelihood {location_to_lik[best_node]}.")
        return best_node, location_to_lik[best_node]

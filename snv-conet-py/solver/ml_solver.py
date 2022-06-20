from typing import Tuple, Set, Optional

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
        self.calculator = LikelihoodCalculator(cell_data=cd, model=model, ctxt=context)

    def insert_snv_events(self) -> generator.SNVModel:
        calculator_refresh = self.calculator.get_likelihood_refresher()
        log_likelihood = calculator_refresh(None)
        logger.info(f"Starting SNV events inference, log-likelihood {log_likelihood}.")
        for snv in self.context.get_snv_event_candidates():
            snv_added = True
            while snv_added:
                snv_added = False
                log_likelihood_snv = calculator_refresh(snv)
                best_node_for_snv, log_lik_with_snv = self.__choose_best_svn_location(snv, calculator_refresh)
                if best_node_for_snv is not None and log_lik_with_snv > log_likelihood_snv:
                    logger.info(f"Adding SNV {snv} to the tree.")
                    self.model.event_tree.node_to_snvs.setdefault(best_node_for_snv, set()).add(snv)
                    calculator_refresh(snv)
                    snv_added = True

        log_likelihood = calculator_refresh(None)
        logger.info(f"Finished SNV events inference. Log-likelihood {log_likelihood}.")
        return self.model

    def __choose_best_svn_location(self, snv: SNVEvent, calculator_refresh) -> Tuple[Optional[CNEvent], float]:
        logger.info(f"Choosing best location for SNV {snv}...")
        location_to_lik = {}
        logger.info(f"Likelihood withut SNV: {calculator_refresh(snv)}")
        for node in self.__get_locations_for_snv(snv):
            self.model.event_tree.node_to_snvs.setdefault(node, set()).add(snv)
            location_to_lik[node] = calculator_refresh(snv)
            logger.info(f"Likelihood for node {node} is {location_to_lik[node]}...")
            self.model.event_tree.node_to_snvs[node].remove(snv)
            calculator_refresh(snv)

        if len(location_to_lik) == 0:
            return None, 0.0
        best_node = [n for n, l in location_to_lik.items() if l == max(location_to_lik.values())][0]
        logger.info(f"Best location is node {best_node} with log-likelihood {location_to_lik[best_node]}.")
        return best_node, location_to_lik[best_node]

    def __get_locations_for_snv(self, snv: SNVEvent) -> Set[CNEvent]:
        return {n for n in self.model.event_tree.cn_event_tree.nodes if self.model.event_tree.snv_can_be_added_in_node(snv, n) and n != EventTreeRoot}
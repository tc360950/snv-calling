from dataclasses import dataclass
from typing import List

import numpy as np
import scipy.stats

import generator
from generator import CNEvent
from solver.cell_data import CellsData

@dataclass
class Genotype:
    altered: int
    cn: int

@dataclass
class LikelihoodCalculator:
    cell_data: CellsData
    model: generator.SNVModel
    context: generator.SNVGeneratorContext


    def get_likelihood(self) -> float: # TODO pzowolic na liczenie dla jednego cell
        d_log_lik = self.__calculate_d_likelihood()
        b_log_lik = self.__calculate_b_likelihood()


    def __calculate_b_likelihood(self) -> np.ndarray:
        b_lik = np.zeros(self.cell_data.d.shape)
        for bin in range(0, b_lik.shape[1]):
            for node in self.model.event_tree.cn_event_tree.nodes:
                for cell in range(0, b_lik.shape[0]):
                    if self.cell_data.attachment[cell] != node:
                        continue
                    genotypes = self.__get_possible_genotypes(node, bin)
                    for g in genotypes:
                        l = scipy.stats.binom.logpmf(self.cell_data.b[cell, bin], self.cell_data.d[cell, bin], max(self.context.get_b_sampling_error(), g.altered / g.cn))
                        # Todo sum


    def __calculate_d_likelihood(self) -> np.ndarray:
        d_lik = np.zeros(self.cell_data.d.shape)
        for cell in range(0, self.cell_data.d.shape[0]):
            node = self.cell_data.attachment[cell]
            for bin in range(0, d_lik.shape[1]):
                if self.model.node_to_cn_profile[node][bin] == 0:
                    d_lik[cell, bin] = \
                        scipy.stats.nbinom.logpmf(self.cell_data.d[cell, bin], self.context.get_total_number_of_reads(),
                                                  1 - self.context.get_b_sampling_error(), loc=0)
                else:
                    d_lik[cell, bin] = \
                        scipy.stats.nbinom.logpmf(self.cell_data.d[cell, bin], self.context.get_total_number_of_reads(), 1 - self.model.node_to_cn_profile[node][bin] * self.context.get_per_allele_coverage(), loc=0)
        return d_lik

    def __get_possible_genotypes(self, node: CNEvent, bin: int) -> List[Genotype]:
        # Todo dla node musisz miec taka macierz ktora mowi czy bin na drodze do node byl altered i czy po altered bylo CN change

        pass
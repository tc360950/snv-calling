from dataclasses import dataclass
from typing import List, Callable

import numpy as np
import scipy.stats
from scipy.special import logsumexp

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

    def get_likelihood(self) -> float:
        d_log_lik: np.ndarray = self.__calculate_d_likelihood()
        b_log_lik: np.ndarray = self.__calculate_b_likelihood()
        for i in range(0, d_log_lik.shape[0]):
            for j in range(0, d_log_lik.shape[1]):
                d_log_lik[i, j] += b_log_lik[i, j]
        return logsumexp(d_log_lik)

    def __calculate_b_likelihood(self) -> np.ndarray:
        def __genotype_log_lik(g: Genotype, b: int, d: int) -> float:
            if g.altered == 0:
                return scipy.stats.binom.logpmf(b, d, self.context.get_b_sampling_error())
            return scipy.stats.binom.logpmf(b, d, g.altered / g.cn)

        b_lik = np.zeros(self.cell_data.d.shape)
        for bin in range(0, b_lik.shape[1]):
            for node in self.model.event_tree.cn_event_tree.nodes:
                genotypes = self.__get_possible_genotypes(node)(bin)
                for cell in range(0, b_lik.shape[0]):
                    if self.cell_data.attachment[cell] != node:
                        continue
                    b_lik[cell, bin] = \
                        logsumexp([-np.log(len(genotypes)) + __genotype_log_lik(g, self.cell_data.b[cell, bin], self.cell_data.d[cell, bin]) for g in genotypes])
        return b_lik

    def __calculate_d_likelihood(self) -> np.ndarray:
        def __get_bin_log_lik(d: int, cn: int) -> float:
            mean = cn * self.context.get_per_allele_coverage()
            if cn == 0:
                mean = self.context.get_d_sampling_error()
            num_failures = int(
                (1.0 - self.context.get_read_success_probablity()) * mean / self.context.get_read_success_probablity())
            return scipy.stats.nbinom.logpmf(d, num_failures, 1 - self.context.get_read_success_probablity())

        d_lik = np.zeros(self.cell_data.d.shape)
        for cell in range(0, self.cell_data.d.shape[0]):
            node = self.cell_data.attachment[cell]
            for bin in range(0, d_lik.shape[1]):
                d_lik[cell, bin] = __get_bin_log_lik(self.cell_data.d[cell, bin], self.model.node_to_cn_profile[node][bin])
        return d_lik

    def __get_possible_genotypes(self, node: CNEvent) -> Callable[[int], List[Genotype]]:
        snvs = self.model.event_tree.node_to_snvs[node]
        bit_map = np.zeros((self.context.get_number_of_bins()))
        self.model.event_tree.mark_bins_with_cn_change_after_alteration(node, bit_map)

        def __get_possible_genotypes_for_bin(bin: int) -> List[Genotype]:
            cn = self.model.node_to_cn_profile[node][bin]
            if bin not in snvs:
                return [Genotype(cn=cn, altered=0)]
            if bit_map[bin] != 1.0:
                return [Genotype(cn=cn, altered=1)]
            return [Genotype(cn=cn, altered=a) for a in range(0, cn)]

        return __get_possible_genotypes_for_bin

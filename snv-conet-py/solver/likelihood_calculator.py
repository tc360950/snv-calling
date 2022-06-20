from dataclasses import dataclass
from typing import List, Callable, Optional, Iterable

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
    ctxt: generator.SNVGeneratorContext

    def get_likelihood_refresher(self):
        b_log_lik = np.zeros(self.cell_data.d.shape)
        d_log_lik = np.zeros(self.cell_data.d.shape)
        self.__calculate_d_likelihood(range(0, self.ctxt.get_number_of_bins()), d_log_lik)
        self.__calculate_b_likelihood(range(0, self.ctxt.get_number_of_bins()), b_log_lik)

        def __calculate(bin: Optional[int]) -> float:
            if bin is None:
                return np.sum(d_log_lik + b_log_lik)
            self.__calculate_b_likelihood([bin], b_log_lik)
            return np.sum(d_log_lik[:, bin] + b_log_lik[:, bin])

        return __calculate

    def __calculate_b_likelihood(self, bins: Iterable[int], b_lik: np.ndarray) -> None:
        def l(g: Genotype, b: int, d: int) -> float:
            if g.altered == 0:
                return scipy.stats.binom.logpmf(b, d, self.ctxt.get_b_sampling_error())
            return scipy.stats.binom.logpmf(b, d, g.altered / g.cn)

        for node in self.model.event_tree.cn_event_tree.nodes:
            genotypes_supplier = self.__get_possible_genotypes_supplier(node)
            for cell in [c for c in range(0, b_lik.shape[0]) if self.cell_data.attachment[c] == node]:
                for bin in bins:
                    genotypes = genotypes_supplier(bin)
                    b_lik[cell, bin] = \
                        logsumexp([-np.log(len(genotypes)) + l(g, self.cell_data.b[cell, bin],
                                                                                self.cell_data.d[cell, bin]) for g in
                                   genotypes])

    def __calculate_d_likelihood(self, bins: Iterable[int], d_lik: np.ndarray) -> None:
        def __get_bin_log_lik(d: int, cn: int, cluster_size: int) -> float:
            mean = cluster_size * (self.ctxt.get_d_sampling_error() if cn == 0 else cn * self.ctxt.get_per_allele_coverage())
            num_failures = int(
                (1.0 - self.ctxt.get_read_success_probablity()) * mean / self.ctxt.get_read_success_probablity())
            return scipy.stats.nbinom.logpmf(d, num_failures, 1 - self.ctxt.get_read_success_probablity())

        for cell in range(0, self.cell_data.d.shape[0]):
            for bin in bins:
                d_lik[cell, bin] = __get_bin_log_lik(d=self.cell_data.d[cell, bin],
                                                     cn=self.model.node_to_cn_profile[self.cell_data.attachment[cell]][bin],
                                                     cluster_size=self.cell_data.cell_cluster_sizes[cell])

    def __get_possible_genotypes_supplier(self, node: CNEvent) -> Callable[[int], List[Genotype]]:
        bit_map = self.model.event_tree.mark_bins_with_cn_change_after_alteration(node, np.zeros((self.ctxt.get_number_of_bins())))
        bins_with_snvs = self.model.event_tree.mark_bins_with_snv(node, np.zeros((self.ctxt.get_number_of_bins())))

        def __get_possible_genotypes_for_bin(bin: int) -> List[Genotype]:
            cn = self.model.node_to_cn_profile[node][bin]
            if bins_with_snvs[bin] == 0.0 or cn == 0:
                return [Genotype(cn=cn, altered=0)]
            if bit_map[bin] != 1.0:
                return [Genotype(cn=cn, altered=1)]
            return [Genotype(cn=cn, altered=a) for a in range(0, cn + 1)]

        return __get_possible_genotypes_for_bin

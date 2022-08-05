from dataclasses import dataclass
from math import exp, log
from typing import Callable, Iterable, List, Optional

import numpy as np
from core.model_data import CellsData, EventTreeWithCounts, Parameters
from core.types import Bin, CNEvent
from scipy.special import logsumexp
from scipy.stats import binom, nbinom


@dataclass
class Genotype:
    altered: int
    cn: int


@dataclass
class LikelihoodCalculator:
    cell_data: CellsData
    tree: EventTreeWithCounts
    p: Parameters

    def get_likelihood_refresher(self):
        b_log_lik = np.zeros(self.cell_data.b.shape)
        d_log_lik = np.zeros(self.cell_data.d.shape)
        self.__calculate_d_likelihood(range(0, self.cell_data.d.shape[1]), d_log_lik)
        self.__calculate_b_likelihood(range(0, self.cell_data.d.shape[1]), b_log_lik)

        def __calculate(bin: Optional[int]) -> float:
            if bin is None:
                return np.sum(d_log_lik + b_log_lik)
            self.__calculate_b_likelihood([bin], b_log_lik)
            return np.sum(d_log_lik[:, bin] + b_log_lik[:, bin])

        return __calculate

    def __calculate_b_likelihood(self, bins: Iterable[int], b_lik: np.ndarray) -> None:
        def lik(g: Genotype, b: int, d: int) -> float:
            return (
                binom.logpmf(b, d, self.p.e)
                if g.altered == 0
                else binom.logpmf(b, d, g.altered / g.cn)
            )

        for node in self.tree.tree.cn_event_tree.nodes:
            genotypes_supplier = self.__get_possible_genotypes_supplier(node)
            for cell in self.cell_data.get_cells_attached_to_node(node):
                for bin in bins:
                    genotypes = genotypes_supplier(bin)
                    b_lik[cell, bin] = logsumexp(
                        [
                            -np.log(len(genotypes))
                            + lik(
                                g,
                                self.cell_data.b[cell, bin],
                                self.cell_data.d[cell, bin],
                            )
                            for g in genotypes
                        ]
                    )

    def __calculate_d_likelihood(self, bins: Iterable[int], d_lik: np.ndarray) -> None:
        def __get_bin_log_lik(d: int, cn: int, cluster_size: int) -> float:
            mean = cluster_size * (self.p.e if cn == 0 else cn * self.p.m)
            num_failures = int(exp(log(1.0 - self.p.q) + log(mean) - log(self.p.q)))
            return (
                0.0
                if num_failures == 0
                else nbinom.logpmf(d, num_failures, 1 - self.p.q)
            )

        for cell in self.cell_data.get_cells():
            for bin in bins:
                d_lik[cell, bin] = __get_bin_log_lik(
                    d=self.cell_data.d[cell, bin],
                    cn=self.tree.node_to_cn_profile[self.cell_data.attachment[cell]][
                        bin
                    ],
                    cluster_size=self.cell_data.cell_cluster_sizes[cell],
                )

    def __get_possible_genotypes_supplier(
        self, node: CNEvent
    ) -> Callable[[Bin], List[Genotype]]:
        bit_map = self.tree.tree.mark_bins_with_cn_change_after_alteration(
            node, np.full((self.cell_data.d.shape[1]), False)
        )
        bins_with_snvs = self.tree.tree.mark_bins_with_snv(
            node, np.full((self.cell_data.d.shape[1]), False)
        )

        def __get_possible_genotypes_for_bin(b: Bin) -> List[Genotype]:
            cn = self.tree.node_to_cn_profile[node][b]
            if not bins_with_snvs[b] or cn == 0:
                return [Genotype(cn=cn, altered=0)]
            if not bit_map[b]:
                return [Genotype(cn=cn, altered=1)]
            return [Genotype(cn=cn, altered=a) for a in range(0, cn + 1)]

        return __get_possible_genotypes_for_bin

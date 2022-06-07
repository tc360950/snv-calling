from abc import ABC
import random
from typing import List
import numpy as np

from generator.gen_utils import SNVEvent, CNEvent, sample_conditionally_with_replacement



class SNVGeneratorContext(ABC):
    def get_epsilon(self) -> float:
        return 0.05

    def get_total_number_of_reads(self) -> int:
        return 100000

    def get_per_allele_coverage(self) -> float:
        return 0.05

    def get_neutral_cn(self) -> int:
        return 2

    def sample_cn_change(self, event: CNEvent, parent_cn_profile: np.ndarray, overlap_bit_map: np.ndarray) -> int:
        """
                Bit map ma 1 jezeli dany bin pojawia sie w cn event w potomku
                :param event:
                :param parent_cn_profile:
                :param overlap_bit_map:
                :return:
                """

        if np.sum(overlap_bit_map) == 0.0: # No bin is present in child nodes- we can do full deletion
            return sample_conditionally_with_replacement(1, lambda: random.randint(-int(np.min(parent_cn_profile)), 2), lambda x: x != 0)[0]

        min_cn_in_bins_with_overlap = np.min(parent_cn_profile[overlap_bit_map == 1.0])
        max_possible_deletion = -min(2, int(min_cn_in_bins_with_overlap) - 1)
        return sample_conditionally_with_replacement(1, lambda: random.randint(max_possible_deletion, 2), lambda x: x != 0)[0]


    def get_number_of_bins(self) -> int:
        return 1000

    def get_cn_event_candidates(self) -> List[CNEvent]:
        return [(a, b) for a in range(0, 1000) for b in range(0, 1000) if a < b and b - a < 50]

    def get_snv_event_candidates(self) -> List[SNVEvent]:
        return [i for i in range(0, 1000)]

    def sample_number_of_snvs_for_edge(self, num_available_snvs: int) -> int:
        return min(num_available_snvs, np.random.poisson(lam=1.0, size=None))

    def get_number_of_alterations(self, cn_before: int, cn_after: int) -> int:
        if cn_after == 0:
            return 0
        if cn_after == cn_before:
            return 1
        elif cn_after > cn_before:
            return 1 if random.randint(0, cn_before) > 0 else cn_after - cn_before
        else:
            return 1 if random.randint(0, cn_before) < cn_after else 0

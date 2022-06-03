from abc import ABC
import random
from typing import List

from networkx.utils.misc import np

from utils import SNVEvent, CNEvent


class SNVGeneratorContext(ABC):

    def get_neutral_cn(self) -> int:
        return 2

    def sample_cn_change(self, event: CNEvent, parent_cn_profile: np.ndarray, overlap_bit_map: np.ndarray) -> int:
        if np.sum(overlap_bit_map) == 0.0:
            return random.randint(-int(np.min(parent_cn_profile)), 2)
        """
        Bit map ma 1 jezeli dany bin pojawia sie w cn event w potomku
        :param event:
        :param parent_cn_profile:
        :param overlap_bit_map:
        :return:
        """
        return 1

    def get_number_of_bins(self) -> int:
        return 1000

    def get_cn_event_candidates(self) -> List[CNEvent]:
        return [(a, b) for a in range(0, 1000) for b in range(0, 1000) if a < b and b - a < 50]

    def get_snv_event_candidates(self) -> List[SNVEvent]:
        return [i for i in range(0, 1000)]

    def sample_number_of_snvs_for_edge(self, num_available_snvs: int) -> int:
        return min(num_available_snvs, np.random.poisson(lam=1.0, size=None))

    def get_number_of_alterations(self, cn_before: int, cn_after: int) -> int:
        if cn_after == cn_before:
            return 1
        elif cn_after > cn_before:
            return 1 if random.randint(0, cn_before) > 0 else cn_after - cn_before
        else:
            return 1 if random.randint(0, cn_before) < cn_after else 0

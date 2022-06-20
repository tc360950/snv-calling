from typing import List, Dict

from generator import EventTree
from .impl.absolute_difference import AbsoluteDifference
from .impl.ancestry_recall import AncestryRecall
from .impl.node_snv_calling import NodeSNVCalling
from .impl.snv_calling import SNVCalling
from stats.statistic import Statistic


class StatisticsCalculator:
    statistics: List[Statistic] = [SNVCalling(), NodeSNVCalling(), AncestryRecall(), AbsoluteDifference()]

    @classmethod
    def calculate_stats(cls, inferred_tree: EventTree, real_tree: EventTree) -> Dict[str, float]:
        all_scores = {}
        for stat in cls.statistics:
            all_scores = {**all_scores, **stat.get_scores(real_tree, inferred_tree)}
        return all_scores

from typing import Dict

from generator import EventTree
from stats.statistic import Statistic


class AbsoluteDifference(Statistic):
    def get_scores(self, real_tree: EventTree, inferred_tree: EventTree) -> Dict[str, float]:
        score = 0
        for n in real_tree.cn_event_tree.nodes:
            score += len(
                real_tree.node_to_snvs.get(n, set()).symmetric_difference(inferred_tree.node_to_snvs.get(n, set())))
        return {"AbsoluteDifference": score}

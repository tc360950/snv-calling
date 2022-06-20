from typing import Dict

from generator import EventTree
from stats.statistic import Statistic

# TODO
class ClusteringPrecision(Statistic):
    def get_scores(self, real_tree: EventTree, inferred_tree: EventTree) -> Dict[str, float]:
        real_node_snvs = {(n, snv) for n in real_tree.cn_event_tree.nodes for snv in
                          real_tree.node_to_snvs.get(n, set())}
        inferred_node_snvs = {(n, snv) for n in inferred_tree.cn_event_tree.nodes for snv in
                              inferred_tree.node_to_snvs.get(n, set())}
        return {"dsa", 0.0}
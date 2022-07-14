from typing import Dict

from core.event_tree import EventTree
from stats.impl.utils import get_sets_precision, get_sets_sensitivity
from stats.statistic import Statistic


class NodeSNVCalling(Statistic):
    def get_scores(self, real_tree: EventTree, inferred_tree: EventTree) -> Dict[str, float]:
        real_node_snvs = {(n, snv) for n in real_tree.cn_event_tree.nodes for snv in
                          real_tree.node_to_snvs.get(n, set())}
        inferred_node_snvs = {(n, snv) for n in inferred_tree.cn_event_tree.nodes for snv in
                              inferred_tree.node_to_snvs.get(n, set())}

        return {
            "NodeSNVCallingPrecision": get_sets_precision(real_node_snvs, inferred_node_snvs),
            "NodeSNVCallingSensitivity": get_sets_sensitivity(real_node_snvs, inferred_node_snvs)
        }

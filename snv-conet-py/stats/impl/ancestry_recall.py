from typing import Dict, Set

import networkx as nx

from generator import EventTree, CNEvent
from stats.statistic import Statistic


class AncestryRecall(Statistic):
    def get_scores(self, real_tree: EventTree, inferred_tree: EventTree) -> Dict[str, float]:
        real_node_snvs = {(n, snv) for n in real_tree.cn_event_tree.nodes for snv in
                          real_tree.node_to_snvs.get(n, set())}
        inferred_node_snvs = {(n, snv) for n in inferred_tree.cn_event_tree.nodes for snv in
                              inferred_tree.node_to_snvs.get(n, set())}

        inferred_snv_to_nodes = dict([(snv[1], set()) for snv in inferred_node_snvs])
        for n, snv in inferred_node_snvs:
            inferred_snv_to_nodes[snv].add(n)

        real_ancestry_pairs = 0
        found_ancestry_pairs = 0
        for n, snv in real_node_snvs:
            for n2, snv2 in real_node_snvs:
                if n == n2 and snv == snv2:
                    continue
                if n in nx.ancestors(real_tree.cn_event_tree, n2):
                    real_ancestry_pairs += 1
                    if self.__has_ancestor_descendant_pair(inferred_snv_to_nodes.get(snv, set()),
                                                           inferred_snv_to_nodes.get(snv2, set()),
                                                           inferred_tree.cn_event_tree):
                        found_ancestry_pairs += 1
                elif n2 in nx.ancestors(real_tree.cn_event_tree, n):
                    real_ancestry_pairs += 1
                    if self.__has_ancestor_descendant_pair(inferred_snv_to_nodes.get(snv2, set()),
                                                           inferred_snv_to_nodes.get(snv, set()),
                                                           inferred_tree.cn_event_tree):
                        found_ancestry_pairs += 1
        return {
            "AncestryRecall": 0.0 if real_ancestry_pairs == 0 else found_ancestry_pairs / real_ancestry_pairs
        }

    def __has_ancestor_descendant_pair(self, nodes1: Set[CNEvent], nodes2: Set[CNEvent], tree: nx.DiGraph) -> bool:
        for n1 in nodes1:
            for n2 in nodes2:
                if n1 in nx.ancestors(tree, n2):
                    return True
        return False

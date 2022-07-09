from dataclasses import dataclass
from typing import Dict, Set, List

import networkx as nx
import numpy as np

from core.types import CNEvent, SNVEvent, EventTreeRoot
from core.utils import apply_to_nodes_in_order


@dataclass
class EventTree:
    cn_event_tree: nx.DiGraph
    node_to_snvs: Dict[CNEvent, Set[SNVEvent]]

    def create_copy_without_snvs(self) -> 'EventTree':
        return EventTree(cn_event_tree=self.cn_event_tree.copy(), node_to_snvs={})

    def get_parent(self, node: CNEvent) -> CNEvent:
        return next(self.cn_event_tree.predecessors(node))

    def snv_can_be_added_in_node(self, snv: SNVEvent, node: CNEvent) -> bool:
        snvs = []
        apply_to_nodes_in_order(node, self.cn_event_tree, lambda n: snvs.extend(self.node_to_snvs.get(n, set())))
        for n in self.__get_path_from_root(node):
            snvs.extend(self.node_to_snvs.get(n, set()))

        return snv not in snvs

    def mark_bins_in_descendant_events(self, node: CNEvent, bit_map: np.ndarray) -> np.ndarray:
        def __mark_overlapping_bins(event: CNEvent) -> None:
            bit_map[[bin for bin in range(event[0], event[1]) if node[0] <= bin < node[1]]] = 1.0

        for descendant in nx.descendants(self.cn_event_tree, node):
            if descendant != node:
                __mark_overlapping_bins(descendant)
        return bit_map

    def mark_bins_with_snv(self, node, bit_map: np.ndarray) -> np.ndarray:
        for n in self.__get_path_from_root(node):
            bit_map[list(self.node_to_snvs.get(n, set()))] = 1.0
        return bit_map

    def mark_bins_with_cn_change_after_alteration(self, node: CNEvent, bit_map: np.ndarray) -> np.ndarray:
        def __mark(n: CNEvent) -> None:
            bit_map[list(self.node_to_snvs.get(n, set()))] = -1.0
            bit_map[[b for b in range(n[0], n[1]) if bit_map[b] == -1.0]] = 1.0

        for n in self.__get_path_from_root(node):
            __mark(n)
        bit_map[bit_map != 1.0] = 0.0
        return bit_map

    def __get_path_from_root(self, node: CNEvent) -> List[CNEvent]:
        if node == EventTreeRoot:
            return []
        return next(nx.all_simple_paths(source=EventTreeRoot, target=node, G=self.cn_event_tree))


from dataclasses import dataclass
from typing import List, Set

import networkx as nx
import numpy as np
from generator.core.types import Bin, CNEvent, EventTreeRoot, SNVEvent, Node


@dataclass
class EventTree:
    event_tree: nx.DiGraph

    def __copy__(self):
        return EventTree(
            event_tree=self.event_tree.copy()
        )

    def get_breakpoint_loci(self) -> List[Bin]:
        return sorted(
            list(
                set(
                    [x[0] for x in self.event_tree.nodes if isinstance(x, tuple)]
                    + [x[1] for x in self.event_tree.nodes if isinstance(x, tuple)]
                )
            )
        )

    def get_first_cn_parent(self, node: Node) -> CNEvent:
        if isinstance(node, CNEvent):
            return node
        return self.get_first_cn_parent(self.get_node_parent(node))

    def get_node_parent(self, node: Node) -> Node:
        return next(self.event_tree.predecessors(node))

    def get_snvs(self) -> list[SNVEvent]:
        return [n for n in self.event_tree.nodes if isinstance(n, int)]

    def mark_bins_in_descendant_events(
            self, node: Node, bin_bit_map: np.ndarray
    ) -> np.ndarray:
        """
        Set @bin_bit_map[bin] to True if CN in the bin is changed in the subtree of node.
        Node is excluded.
        """

        def __mark_overlapping_bins(event: CNEvent) -> None:
            bin_bit_map[
                [b for b in range(event[0], event[1]) if node[0] <= b < node[1]]
            ] = True

        for descendant in nx.descendants(self.event_tree, node):
            if descendant != node and isinstance(descendant, tuple):
                __mark_overlapping_bins(descendant)
        return bin_bit_map

    def mark_bins_with_snv(self, node, bin_bit_map: np.ndarray) -> np.ndarray:
        bin_bit_map[list(self.__gather_snvs_on_path_to_root(node))] = True
        return bin_bit_map

    def mark_bins_with_cn_change_after_alteration(
            self, node: Node, bin_bit_map: np.ndarray
    ) -> np.ndarray:
        bit_map = np.full(bin_bit_map.shape, 0, dtype=int)
        BIN_IN_SNV = 1
        CN_ALTERATION_AFTER_SNV = 2

        for n in self.__get_path_from_root(node):
            if isinstance(n, int):
                bit_map[n] = BIN_IN_SNV
            else:
                bit_map[
                    [b for b in range(n[0], n[1]) if bit_map[b] == BIN_IN_SNV]
                ] = CN_ALTERATION_AFTER_SNV

        bin_bit_map[bit_map == CN_ALTERATION_AFTER_SNV] = True
        return bit_map

    def __gather_snvs_in_subtree(self, node: Node) -> Set[SNVEvent]:
        return {d for d in nx.descendants(self.event_tree, node) if isinstance(d, int)}

    def gather_snvs_on_path_to_root(self, node: Node) -> Set[SNVEvent]:
        snvs = []
        for n in self.__get_path_from_root(node):
            if isinstance(n, int):
                snvs.append(n)
        return set(snvs)

    def __get_path_from_root(self, node: Node) -> List[Node]:
        return (
            []
            if node == EventTreeRoot
            else next(
                nx.all_simple_paths(
                    source=EventTreeRoot, target=node, G=self.event_tree
                )
            )
        )

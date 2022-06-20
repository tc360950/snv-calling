import random
from dataclasses import dataclass
from typing import Dict, Set, List

import networkx as nx
import numpy as np

from generator.context import SNVGeneratorContext
from generator.tree_sampler import RandomWalkTreeSampler
from generator.gen_utils import CNEvent, SNVEvent, EventTreeRoot, sample_conditionally_without_replacement, \
    apply_to_nodes_in_order


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

    def mark_bins_wth_snv(self, node, bit_map: np.ndarray) -> np.ndarray:
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


@dataclass
class EventTreeGenerator:
    context: SNVGeneratorContext

    def generate(self, tree_size: int) -> EventTree:
        event_tree = RandomWalkTreeSampler.sample_tree(self.__generate_nodes(tree_size))
        return EventTree(event_tree, self.__generate_snvs(event_tree))

    def __generate_nodes(self, tree_size: int) -> List[CNEvent]:
        return [EventTreeRoot] + list(random.sample(self.context.get_cn_event_candidates(), tree_size))

    def __sample_snvs_for_node(self, ancestor_snvs: Set[SNVEvent]) -> Set[SNVEvent]:
        snv_candidates = self.context.get_snv_event_candidates()
        return sample_conditionally_without_replacement(
            k=self.context.sample_number_of_snvs_for_edge(len(snv_candidates) - len(ancestor_snvs)),
            sampler=lambda: random.sample(snv_candidates, 1)[0], condition=lambda x: x not in ancestor_snvs)

    def __generate_snvs(self, cn_tree: nx.DiGraph) -> Dict[CNEvent, Set[SNVEvent]]:
        node_to_snvs = {}

        def __generate_snv_events_for_node(node: CNEvent, ancestor_snvs: Set[SNVEvent]) -> None:
            node_to_snvs[node] = self.__sample_snvs_for_node(ancestor_snvs)
            ancestors_snvs = ancestor_snvs.union(node_to_snvs[node])
            for child in cn_tree.successors(node):
                __generate_snv_events_for_node(child, ancestors_snvs)

        __generate_snv_events_for_node(EventTreeRoot, set())
        return node_to_snvs


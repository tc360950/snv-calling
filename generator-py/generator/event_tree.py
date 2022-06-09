import random
from dataclasses import dataclass
from typing import Dict, Set, List

import networkx as nx
import numpy as np

from generator.context import SNVGeneratorContext
from generator.tree_sampler import RandomWalkTreeSampler
from generator.gen_utils import CNEvent, SNVEvent, EventTreeRoot, sample_conditionally_without_replacement


@dataclass
class EventTree:
    cn_event_tree: nx.DiGraph
    node_to_snvs: Dict[CNEvent, Set[SNVEvent]]

    def get_parent(self, node: CNEvent) -> CNEvent:
        return next(self.cn_event_tree.predecessors(node))

    def mark_overlapping_bins_from_descendant(self, node: CNEvent, bit_map: np.ndarray) -> None:
        for descendant in nx.descendants(self.cn_event_tree, node):
            if descendant != node:
                self.__mark_overlapping_bins(descendant, node, bit_map)

    def mark_bins_with_cn_change_after_alteration(self, node: CNEvent, bit_map: np.ndarray) -> None:
        def __mark(n: CNEvent) -> None:
            if n == EventTreeRoot:
                return
            bit_map[list(self.node_to_snvs[n])] = -1.0
            for bin in range(n[0], n[1]):
                if bit_map[bin] == -1.0:
                    bit_map[bin] = 1.0
        for n in self.__get_path_from_root(node):
            __mark(n)

        for bin in range(0, bit_map.shape[0]):
            if bit_map[bin] != 1.0:
                bit_map[bin] = 0.0


    def __get_path_from_root(self, node: CNEvent) -> List[CNEvent]:
        result = [node]
        while node != EventTreeRoot:
            result.insert(0, node)
        result.insert(0, EventTreeRoot)
        return result

    def __mark_overlapping_bins(self, event: CNEvent, reference_event: CNEvent, bit_map: np.ndarray) -> None:
        for bin in range(event[0], event[1]):
            if reference_event[0] <= bin < reference_event[1]:
                bit_map[bin] = 1.0


@dataclass
class EventTreeGenerator:
    context: SNVGeneratorContext

    def generate(self, tree_size: int) -> EventTree:
        cn_nodes = [EventTreeRoot] + list(random.sample(self.context.get_cn_event_candidates(), tree_size))
        event_tree = RandomWalkTreeSampler.sample_tree(cn_nodes)
        node_to_snvs = {}
        self.__generate_snv_events_for_node(EventTreeRoot, event_tree, set(), node_to_snvs)
        return EventTree(event_tree, node_to_snvs)

    def __generate_snv_events_for_node(self, tree_node: CNEvent, cn_tree: nx.DiGraph, ancestors_snvs: Set[SNVEvent],
                                       node_to_snvs: Dict[CNEvent, Set[SNVEvent]]) -> None:
        snv_candidates = self.context.get_snv_event_candidates()
        events = sample_conditionally_without_replacement(
            k=self.context.sample_number_of_snvs_for_edge(len(snv_candidates) - len(ancestors_snvs)),
            sampler=lambda: random.sample(snv_candidates, 1)[0], condition=lambda x: x not in ancestors_snvs)
        node_to_snvs[tree_node] = events
        ancestors_snvs = ancestors_snvs.union(events)
        for child in cn_tree.successors(tree_node):
            self.__generate_snv_events_for_node(child, cn_tree, ancestors_snvs, node_to_snvs)

import random
from dataclasses import dataclass
from typing import Dict, Set, List

import networkx as nx

from core.event_tree import EventTree
from core.types import CNEvent, SNVEvent, EventTreeRoot
from generator.context import SNVGeneratorContext
from generator.tree_sampler import RandomWalkTreeSampler
from generator.gen_utils import sample_conditionally_without_replacement

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
            if node != EventTreeRoot:
                node_to_snvs[node] = self.__sample_snvs_for_node(ancestor_snvs)
                ancestor_snvs = ancestor_snvs.union(node_to_snvs[node])
            for child in cn_tree.successors(node):
                __generate_snv_events_for_node(child, ancestor_snvs)

        __generate_snv_events_for_node(EventTreeRoot, set())
        return node_to_snvs


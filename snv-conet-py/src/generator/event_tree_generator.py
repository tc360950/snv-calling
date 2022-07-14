import random
from dataclasses import dataclass
from typing import Dict, Set, List, TypeVar

import networkx as nx

from core.event_tree import EventTree
from core.types import CNEvent, SNVEvent, EventTreeRoot
from generator.context import SNVGeneratorContext
from generator.gen_utils import sample_conditionally_without_replacement

T = TypeVar('T')


class RandomWalkTreeSampler:
    @staticmethod
    def sample_tree(nodes: List[T]) -> nx.DiGraph:
        """
            Sample random uniform directed rooted tree built from (all) nodes in @nodes.
            Always nodes[0] is assumed to be the root.
        """
        tree = nx.DiGraph()
        tree.add_node(nodes[0])
        current_node = nodes[0]
        while len(nodes) > tree.size() + 1:  # nx.DiGraph.size evaluates to the number of graph edges
            new_node = random.sample(range(0, len(nodes)), 1)[0]
            if nodes[new_node] not in set(tree.nodes):
                tree.add_edge(current_node, nodes[new_node])
            current_node = nodes[new_node]
        return tree


@dataclass
class EventTreeGenerator:
    context: SNVGeneratorContext

    def generate(self, tree_size: int) -> EventTree:
        event_tree = RandomWalkTreeSampler.sample_tree(self.__generate_cn_event_nodes(tree_size))
        return EventTree(event_tree, self.__generate_snvs(event_tree))

    def __generate_cn_event_nodes(self, tree_size: int) -> List[CNEvent]:
        return [EventTreeRoot] + list(random.sample(self.context.get_cn_event_candidates(), tree_size))

    def __generate_snvs(self, cn_tree: nx.DiGraph) -> Dict[CNEvent, Set[SNVEvent]]:
        node_to_snvs: Dict[CNEvent, Set[SNVEvent]] = {}

        def generate_snv_events_for_node(node: CNEvent, ancestor_snvs: Set[SNVEvent]) -> None:
            if node != EventTreeRoot:
                node_to_snvs[node] = self.__sample_snvs_for_node(ancestor_snvs)
                ancestor_snvs = ancestor_snvs.union(node_to_snvs[node])
            for child in cn_tree.successors(node):
                generate_snv_events_for_node(child, ancestor_snvs)

        generate_snv_events_for_node(EventTreeRoot, set())
        return node_to_snvs

    def __sample_snvs_for_node(self, ancestor_snvs: Set[SNVEvent]) -> Set[SNVEvent]:
        snv_candidates = self.context.get_snv_event_candidates()
        return sample_conditionally_without_replacement(
            k=self.context.sample_number_of_snvs_for_edge(len(snv_candidates) - len(ancestor_snvs)),
            sampler=lambda: random.sample(snv_candidates, 1)[0],
            condition=lambda x: x not in ancestor_snvs)

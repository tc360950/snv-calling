import random
from dataclasses import dataclass
from typing import List, TypeVar

import networkx as nx
from generator.core.event_tree import EventTree
from generator.core.types import EventTreeRoot, Node
from generator.generator.context import SNVGeneratorContext

T = TypeVar("T")


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
        while (
            len(nodes) > tree.size() + 1
        ):  # nx.DiGraph.size evaluates to the number of graph edges
            new_node = random.sample(range(0, len(nodes)), 1)[0]
            if nodes[new_node] not in set(tree.nodes):
                tree.add_edge(current_node, nodes[new_node])
            current_node = nodes[new_node]
        return tree


@dataclass
class EventTreeGenerator:
    context: SNVGeneratorContext

    def generate(self, tree_size: int) -> EventTree:
        event_tree = RandomWalkTreeSampler.sample_tree(
            self.__generate_nodes(tree_size)
        )
        return EventTree(event_tree)

    def __generate_nodes(self, tree_size: int) -> List[Node]:
        snvs_count = 0

        for _ in range(1, tree_size):
            if random.random() <= self.context.snv_event_prob():
                snvs_count += 1
        cn_nodes = list(
            random.sample(self.context.get_cn_event_candidates(), tree_size - snvs_count - 1)
        )
        snv_nodes = list(
            random.sample(self.context.get_snv_event_candidates(), snvs_count)
        )
        nodes = cn_nodes + snv_nodes
        random.shuffle(nodes)
        return [EventTreeRoot] + nodes


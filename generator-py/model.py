import random
from typing import Generator, Any, Set

import numpy as np

from context import SNVGeneratorContext
from event_tree import EventTreeGenerator, EventTree
from tree_sampler import RandomWalkTreeSampler
from utils import CNEvent, SNVEvent, EventTreeRoot, sample_conditionally, cn_events_overlap


class SNVModel:
    context: SNVGeneratorContext

    def __init__(self, context: SNVGeneratorContext):
        self.context = context
        self.event_tree: EventTree = None
        self.node_to_cn_profile = {EventTreeRoot: np.full((context.get_number_of_bins()), fill_value=context.get_neutral_cn())}
        self.node_to_altered_counts_profile = {EventTreeRoot: np.zeros((context.get_number_of_bins()))}

    def generate_structure(self, tree_size: int) -> None:
        self.event_tree = EventTreeGenerator(self.context).generate(tree_size)
        for child in self.event_tree.cn_event_tree.successors(EventTreeRoot):
            self.__generate_cn_profile(child)
            self.__generate_altered_counts_profile(child)

    def __generate_cn_profile(self, tree_node: CNEvent) -> None:
        parent_cn_profile = self.node_to_cn_profile[self.event_tree.get_parent(tree_node)].copy()
        overlap_bitmap = np.zeros((self.context.get_number_of_bins()))
        cn_change = self.context.sample_cn_change(tree_node, parent_cn_profile, self.event_tree.mark_overlapping_bins_from_descendant(tree_node, overlap_bitmap))
        for bin in range(tree_node[0], tree_node[1]):
            parent_cn_profile[bin] += cn_change
        self.node_to_cn_profile[tree_node] = parent_cn_profile
        for child in self.event_tree.cn_event_tree.successors(tree_node):
            self.__generate_cn_profile(child)

    def __generate_altered_counts_profile(self, node: CNEvent) -> None:
        parent_alteration_counts = self.node_to_altered_counts_profile[self.event_tree.get_parent(node)]
        node_cn_profile = self.node_to_cn_profile[node]

        for snv in self.event_tree.node_to_snvs[node]:
            assert parent_alteration_counts[snv] == 0
            parent_alteration_counts[snv] = self.context.get_number_of_alterations(self.node_to_cn_profile[self.event_tree.get_parent(node)][snv], node_cn_profile[snv])

        self.node_to_altered_counts_profile[node] = parent_alteration_counts
        for child in self.event_tree.cn_event_tree.successors(node):
            self.__generate_altered_counts_profile(child)





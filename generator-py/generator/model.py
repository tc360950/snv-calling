from dataclasses import dataclass
from typing import Set, Generator, Any, Dict

import networkx as nx
import numpy as np
import random

from generator.context import SNVGeneratorContext
from generator.event_tree import EventTreeGenerator, EventTree
from generator.gen_utils import CNEvent, EventTreeRoot


@dataclass
class CellData:
    d: np.ndarray
    b: np.ndarray
    attachment: CNEvent


class SNVModel:
    context: SNVGeneratorContext

    def __init__(self, context: SNVGeneratorContext):
        self.context = context
        self.event_tree: EventTree = None
        self.node_to_cn_profile = {
            EventTreeRoot: np.full((context.get_number_of_bins()), fill_value=context.get_neutral_cn())}
        self.node_to_altered_counts_profile = {EventTreeRoot: np.zeros((context.get_number_of_bins()))}

    def generate_structure(self, tree_size: int) -> None:
        self.event_tree = EventTreeGenerator(self.context).generate(tree_size)
        for child in self.event_tree.cn_event_tree.successors(EventTreeRoot):
            self.__generate_cn_profile(child)
            self.__generate_altered_counts_profile(child)

    def generate_cell(self) -> Generator[CellData, CellData, Any]:
        if self.event_tree is None:
            raise RuntimeError(f"Can;'t generate cells withput structure genmeration. Call generate_structure first!")

        probs = self.__create_node_attachment_probabilities()

        while True:
            yield self.__sample_cell_data(random.choices(list(probs), k=1, weights=list(probs.values()))[0])

    def __generate_cn_profile(self, tree_node: CNEvent) -> None:
        parent_cn_profile = self.node_to_cn_profile[self.event_tree.get_parent(tree_node)].copy()
        overlap_bitmap = np.zeros((self.context.get_number_of_bins()))
        self.event_tree.mark_overlapping_bins_from_descendant(tree_node, overlap_bitmap)
        cn_change = self.context.sample_cn_change(tree_node, parent_cn_profile, overlap_bitmap)
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
            parent_alteration_counts[snv] = self.context.get_number_of_alterations(
                self.node_to_cn_profile[self.event_tree.get_parent(node)][snv], node_cn_profile[snv])

        self.node_to_altered_counts_profile[node] = parent_alteration_counts
        for child in self.event_tree.cn_event_tree.successors(node):
            self.__generate_altered_counts_profile(child)

    def __sample_cell_data(self, node: CNEvent) -> CellData:
        # TODO co z CN = 0 dla b i d?
        d = np.zeros(self.context.get_number_of_bins())
        b = np.zeros(self.context.get_number_of_bins())


        for bin in range(0, d.shape[0]):
            if self.node_to_cn_profile[node][bin] == 0:
                d[bin] = np.random.binomial(self.context.get_total_number_of_reads(), self.context.get_epsilon())
            else:
                d[bin] = np.random.negative_binomial(self.context.get_total_number_of_reads(), self.node_to_cn_profile[node][bin] * self.context.get_per_allele_coverage())

            if self.node_to_altered_counts_profile[node][bin] == 0:
                b[bin] = np.random.binomial(d[bin], self.context.get_epsilon())
            else:
                b[bin] = np.random.binomial(d[bin], self.node_to_altered_counts_profile[node][bin] / self.node_to_cn_profile[node][bin])

        return CellData(d=d, b =b, attachment=node)

    def __create_node_attachment_probabilities(self) -> Dict[CNEvent, float]:
        """
        Probability of cell being attached to given node is proportional to  e^{0.1 * depth} where depth is node's depth in the tree.
        """
        depths = nx.shortest_path_length(self.event_tree.cn_event_tree, source=EventTreeRoot)
        depths = dict([(key, np.exp(0.1 * val)) for key, val in depths.items()])
        sum_depths = sum(list(depths.values()))
        for key, value in depths.items():
            depths[key] = np.exp(0.1 * depths[key]) / sum_depths
        return depths


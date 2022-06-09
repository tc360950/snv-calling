from dataclasses import dataclass
from typing import Set, Generator, Any, Dict

import networkx as nx
import numpy as np
import random

from generator.context import SNVGeneratorContext
from generator.event_tree import EventTreeGenerator, EventTree
from generator.gen_utils import CNEvent, EventTreeRoot, apply_to_nodes_in_order


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
        apply_to_nodes_in_order(EventTreeRoot, self.event_tree.cn_event_tree, lambda n: self.__generate_node_cn_profile(n))
        apply_to_nodes_in_order(EventTreeRoot, self.event_tree.cn_event_tree, lambda n: self.__generate_node_altered_counts_profile(n))

    def generate_cell(self) -> Generator[CellData, CellData, Any]:
        if self.event_tree is None:
            raise RuntimeError(f"Can;'t generate cells without structure generation. Call generate_structure first!")

        probs = self.__create_node_attachment_probabilities()

        while True:
            yield self.__sample_cell_data(random.choices(list(probs), k=1, weights=list(probs.values()))[0])

    def __generate_node_cn_profile(self, node: CNEvent) -> None:
        if node == EventTreeRoot:
            return
        parent_cn_profile = self.node_to_cn_profile[self.event_tree.get_parent(node)].copy()
        overlap_bitmap = np.zeros((self.context.get_number_of_bins()))
        self.event_tree.mark_overlapping_bins_from_descendant(node, overlap_bitmap)
        cn_change = self.context.sample_cn_change(node, parent_cn_profile, overlap_bitmap)
        for bin in range(node[0], node[1]):
            parent_cn_profile[bin] += cn_change
        self.node_to_cn_profile[node] = parent_cn_profile


    def __generate_node_altered_counts_profile(self, node: CNEvent) -> None:
        if node == EventTreeRoot:
            return
        parent_alteration_counts = self.node_to_altered_counts_profile[self.event_tree.get_parent(node)]
        node_cn_profile = self.node_to_cn_profile[node]

        for snv in self.event_tree.node_to_snvs[node]:
            parent_alteration_counts[snv] = self.context.get_number_of_alterations(
                self.node_to_cn_profile[self.event_tree.get_parent(node)][snv], node_cn_profile[snv])

        self.node_to_altered_counts_profile[node] = parent_alteration_counts

    def __sample_read_in_bin(self, cn: int) -> int:
        mean = cn * self.context.get_per_allele_coverage()
        if cn == 0:
            mean = self.context.get_d_sampling_error()
        num_failures = int(
            (1 - self.context.get_read_success_probablity()) * mean / self.context.get_read_success_probablity())
        return np.random.negative_binomial(num_failures, self.context.get_read_success_probablity())

    def __sample_altered_reads_in_bin(self, altered_counts: int, total_reads: int, cn: int) -> int:
        if altered_counts == 0:
            return np.random.binomial(total_reads, self.context.get_b_sampling_error())
        return np.random.binomial(total_reads, altered_counts / cn)

    def __sample_cell_data(self, node: CNEvent) -> CellData:
        d = np.zeros(self.context.get_number_of_bins())
        b = np.zeros(self.context.get_number_of_bins())

        for bin in range(0, d.shape[0]):
            d[bin] = self.__sample_read_in_bin(self.node_to_cn_profile[node][bin])
            b[bin] = self.__sample_altered_reads_in_bin(self.node_to_altered_counts_profile[node][bin], d[bin], self.node_to_cn_profile[node][bin])

        return CellData(d=d, b=b, attachment=node)

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

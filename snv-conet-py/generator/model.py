from dataclasses import dataclass
from typing import Generator, Any, Dict

import networkx as nx
import numpy as np
import random

from core.types import CNEvent, EventTreeRoot
from core.utils import apply_to_nodes_in_order
from generator.cn_sampler import CNSampler
from generator.context import SNVGeneratorContext
from generator.corrected_counts_generator import CountsGenerator
from generator.event_tree_generator import EventTreeGenerator, EventTree


@dataclass
class CellData:
    d: np.ndarray
    b: np.ndarray
    attachment: CNEvent
    corrected_counts: np.ndarray


class SNVModel:
    ctxt: SNVGeneratorContext
    counts_gen: CountsGenerator

    def __init__(self, ctxt: SNVGeneratorContext):
        self.ctxt = ctxt
        self.event_tree: EventTree = None
        self.node_to_cn_profile = {
            EventTreeRoot: np.full((ctxt.get_number_of_bins()), fill_value=ctxt.get_neutral_cn())}
        self.node_to_altered_counts_profile = {EventTreeRoot: np.zeros((ctxt.get_number_of_bins()))}
        self.counts_gen = CountsGenerator(CNSampler.create_default_sampler())

    def generate_structure(self, tree_size: int) -> None:
        self.event_tree = EventTreeGenerator(self.ctxt).generate(tree_size)
        apply_to_nodes_in_order(EventTreeRoot, self.event_tree.cn_event_tree,
                                lambda n: self.__generate_node_cn_profile(n) if n != EventTreeRoot else None)
        apply_to_nodes_in_order(EventTreeRoot, self.event_tree.cn_event_tree,
                                lambda n: self.__generate_node_altered_counts_profile(n) if n != EventTreeRoot else None)

    def generate_cell_in_node(self, node: CNEvent) -> CellData:
        return self.__sample_cell_data(node)

    def generate_cell_attachment(self) -> Generator[CNEvent, CNEvent, Any]:
        if self.event_tree is None:
            raise RuntimeError(f"Can't generate cells without structure generation. Call generate_structure first!")

        probs = self.__create_node_attachment_probabilities()
        while True:
            yield random.choices(list(probs), k=1, weights=list(probs.values()))[0]

    def __generate_node_cn_profile(self, node: CNEvent) -> None:
        self.node_to_cn_profile[node] = self.node_to_cn_profile[self.event_tree.get_parent(node)].copy()
        bitmap = self.event_tree.mark_bins_in_descendant_events(node, np.zeros((self.ctxt.get_number_of_bins())))
        cn_change = self.ctxt.sample_cn_change(node, self.node_to_cn_profile[self.event_tree.get_parent(node)], bitmap)
        for bin in range(node[0], node[1]):
            self.node_to_cn_profile[node][bin] += cn_change

    def __generate_node_altered_counts_profile(self, node: CNEvent) -> None:
        parent_alteration_counts = self.node_to_altered_counts_profile[self.event_tree.get_parent(node)].copy()
        node_cn_profile = self.node_to_cn_profile[node]

        for snv in self.event_tree.node_to_snvs[node]:
            if self.node_to_cn_profile[self.event_tree.get_parent(node)][snv] > 0:
                parent_alteration_counts[snv] = 1 # TODO jakos do context
        for bin in range(node[0], node[1]):
            if parent_alteration_counts[bin] > 0.0:
                parent_alteration_counts[bin] = \
                    self.ctxt.get_number_of_alterations(self.node_to_cn_profile[self.event_tree.get_parent(node)][bin],
                                                        node_cn_profile[bin], parent_alteration_counts[bin])

        self.node_to_altered_counts_profile[node] = parent_alteration_counts

    def __sample_read_in_bin(self, cn: int) -> int:
        mean = cn * self.ctxt.get_per_allele_coverage()
        if cn == 0:
            mean = self.ctxt.get_d_sampling_error()
        num_failures = int(
            (1.0 - self.ctxt.get_read_success_probablity()) * mean / self.ctxt.get_read_success_probablity())
        return np.random.negative_binomial(num_failures, 1 - self.ctxt.get_read_success_probablity())

    def __sample_altered_reads_in_bin(self, altered_counts: int, total_reads: int, cn: int) -> int:
        if total_reads == 0:
            return 0
        if altered_counts == 0 or cn == 0:
            return np.random.binomial(total_reads, self.ctxt.get_b_sampling_error())
        return np.random.binomial(total_reads, altered_counts / cn)

    def __sample_cell_data(self, node: CNEvent) -> CellData:
        d = np.zeros(self.ctxt.get_number_of_bins())
        b = np.zeros(self.ctxt.get_number_of_bins())

        for bin in range(0, d.shape[0]):
            d[bin] = self.__sample_read_in_bin(self.node_to_cn_profile[node][bin])
            b[bin] = self.__sample_altered_reads_in_bin(self.node_to_altered_counts_profile[node][bin], d[bin],
                                                        self.node_to_cn_profile[node][bin])

        return CellData(d=d, b=b, attachment=node, corrected_counts= self.counts_gen.add_noise_to_counts(self.node_to_cn_profile[node]))

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

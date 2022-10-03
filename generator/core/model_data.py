from dataclasses import dataclass
from typing import Dict, List

import numpy as np
from generator.core.event_tree import EventTree
from generator.core.types import Cell, CNEvent


@dataclass
class Parameters:
    m: float  # per allele coverage
    e: float  # sequencing error
    q: float  # read success probablity

    def __str__(self) -> str:
        return (
            f"Per allele coverage: {self.m}, "
            f"sequencing error: {self.e}, "
            f"read success probability: {self.q}."
        )


@dataclass
class EventTreeWithCounts:
    tree: EventTree
    node_to_cn_profile: Dict[CNEvent, np.ndarray]


@dataclass
class CellsData:
    d: np.ndarray
    b: np.ndarray
    attachment: List[CNEvent]
    # Each cell represents a cluster of given size
    cell_cluster_sizes: List[int]

    def get_cells_attached_to_node(self, node: CNEvent) -> List[Cell]:
        return [c for c in range(0, len(self.attachment)) if self.attachment[c] == node]

    def get_cells(self) -> List[Cell]:
        return list(range(0, len(self.attachment)))

    def aggregate(self) -> "CellsData":
        """
        Collapse same node clusters into one big cluster
        """
        collapsed_d = np.zeros((len(set(self.attachment)), self.d.shape[1]))
        collapsed_b = np.zeros((len(set(self.attachment)), self.b.shape[1]))
        new_cluster_sizes = []
        new_attachment = []
        counter = 0
        for node in set(self.attachment):
            same_node_clusters = [
                i for i in range(0, len(self.attachment)) if self.attachment[i] == node
            ]
            new_cluster_sizes.append(
                sum([self.cell_cluster_sizes[i] for i in same_node_clusters])
            )
            new_attachment.append(node)
            for i in same_node_clusters:
                collapsed_d[counter, :] += self.d[i, :]
                collapsed_b[counter, :] += self.b[i, :]
            counter += 1
        return CellsData(
            d=collapsed_d,
            b=collapsed_b,
            attachment=new_attachment,
            cell_cluster_sizes=new_cluster_sizes,
        )


def create_cn_matrix(cells: CellsData, tree: EventTreeWithCounts) -> np.ndarray:
    matrix = np.zeros(cells.d.shape)
    for cell in range(0, cells.d.shape[0]):
        matrix[
            cell,
        ] = tree.node_to_cn_profile[cells.attachment[cell]]
    return matrix.astype(int)

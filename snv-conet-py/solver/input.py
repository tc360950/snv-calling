from dataclasses import dataclass
from typing import Dict

import numpy as np

from core.event_tree import EventTree
from core.types import CNEvent


@dataclass
class EventTreeWithCounts:
    tree: EventTree
    node_to_cn_profile: Dict[CNEvent, np.ndarray]


@dataclass
class CellsData:
    d: np.ndarray
    b: np.ndarray
    attachment: Dict[int, CNEvent]
    # Each cell represents a cluster of given size
    cell_cluster_sizes: Dict[int, int]


def create_cn_matrix(cells: CellsData, tree: EventTreeWithCounts) -> np.ndarray:
    matrix = np.zeros(cells.d.shape)
    for cell in range(0, cells.d.shape[0]):
        matrix[cell,] = tree.node_to_cn_profile[cells.attachment[cell]]
    return matrix.astype(int)

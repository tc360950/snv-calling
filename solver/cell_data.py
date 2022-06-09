from dataclasses import dataclass
from typing import Dict

import numpy as np

import generator.gen_utils


@dataclass
class CellsData:
    d: np.ndarray
    b: np.ndarray
    attachment: Dict[int, generator.gen_utils.CNEvent]
from abc import ABC, abstractmethod
from typing import Dict

from generator import EventTree


class Statistic(ABC):
    @abstractmethod
    def get_scores(self, real_tree: EventTree, inferred_tree: EventTree) -> Dict[str, float]:
        pass
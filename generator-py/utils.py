import logging
import sys
from typing import Tuple, TypeVar, Callable, Set

CNEvent = Tuple[int, int]
SNVEvent = int
EventTreeRoot: CNEvent = (0, 0)

T = TypeVar('T')


def cn_events_overlap(ev1: CNEvent, ev2: CNEvent) -> bool:
    return (ev1[0] <= ev2[0] < ev1[1]) or (ev2[0] <= ev1[0] < ev2[1])


def __sample_conditionally(sampler: Callable[[], T], condition: Callable[[T], bool]) -> T:
    """
        Returns first result of call to @sampler which satisfies @condition
    """
    sample = sampler()
    while not condition(sample):
        sample = sampler()
    return sample


def sample_conditionally(k: int, sampler: Callable[[], T], condition: Callable[[T], bool]) -> Set[T]:
    result = set()
    for _ in range(0, k):
        result.add(__sample_conditionally(sampler, lambda x: x not in result and condition(x)))
    return result


FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s — %(message)s")


def get_logger(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(FORMATTER)
    logger.addHandler(handler)
    logger.propagate = False
    return logger

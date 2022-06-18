from typing import TypeVar, Set

T = TypeVar('T')


def get_sets_precision(real_set: Set[T], inferred_set: Set[T]) -> float:
    return 0.0 if not inferred_set else len(real_set.intersection(inferred_set)) / len(inferred_set)


def get_sets_sensitivity(real_set: Set[T], inferred_set: Set[T]) -> float:
    return 0.0 if not real_set else len(real_set.intersection(inferred_set)) / len(real_set)

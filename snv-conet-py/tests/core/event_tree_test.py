import networkx as nx
import numpy as np
import pytest as pytest

from core.event_tree import EventTree
from core.types import EventTreeRoot


@pytest.fixture
def event_tree() -> EventTree:
    tree = nx.DiGraph()
    tree.add_node(EventTreeRoot)
    tree.add_edge(EventTreeRoot, (0, 10))
    tree.add_edge((0, 10), (5, 15))

    return EventTree(cn_event_tree=tree, node_to_snvs={(0, 10): {12}, (5, 15): {11}})


def test_snv_addition(event_tree: EventTree):
    assert event_tree.snv_can_be_added_in_node(0, (0, 10))
    assert not event_tree.snv_can_be_added_in_node(11, (0, 10))
    assert not event_tree.snv_can_be_added_in_node(12, (5, 15))


def test_snv_bitmaps(event_tree: EventTree):
    bit_map = np.full(20, False)
    event_tree.mark_bins_with_snv((0, 10), bit_map)
    assert bit_map[12]
    for i in range(0, 20):
        assert not bit_map[i] or i == 12

    bit_map = event_tree.mark_bins_with_snv((5, 15), np.full(20, False))
    assert bit_map[12]
    assert bit_map[11]
    for i in range(0, 20):
        assert not bit_map[i] or i in {12, 11}


def test_bins_with_cn_change_bitmap(event_tree: EventTree):
    bit_map = np.full(20, False)
    event_tree.mark_bins_with_cn_change_after_alteration((0, 10), bit_map)
    assert not np.any(bit_map)

    bit_map = event_tree.mark_bins_with_cn_change_after_alteration((5, 15), np.full(20, False))
    assert bit_map[11] and bit_map[12]
    for i in range(0, 20):
        assert not bit_map[i] or i in {12, 11}


def test_bins_in_descendant_events_bitmap(event_tree: EventTree):
    bit_map = event_tree.mark_bins_in_descendant_events((0, 10), np.full(20, False))
    assert np.all(bit_map[5:10])

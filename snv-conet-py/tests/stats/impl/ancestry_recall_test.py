import networkx as nx
import pytest

from core.event_tree import EventTree
from core.types import EventTreeRoot
from stats.impl.ancestry_recall import AncestryRecall


@pytest.fixture
def event_tree() -> EventTree:
    tree = nx.DiGraph()
    tree.add_node(EventTreeRoot)
    tree.add_edge(EventTreeRoot, (0, 10))
    tree.add_edge((0, 10), (5, 15))

    return EventTree(cn_event_tree=tree, node_to_snvs={(0, 10): {0, 12}, (5, 15): {1, 11}})


@pytest.fixture
def event_tree2() -> EventTree:
    tree = nx.DiGraph()
    tree.add_node(EventTreeRoot)
    tree.add_edge(EventTreeRoot, (0, 10))
    tree.add_edge((0, 10), (5, 15))
    tree.add_edge((0, 10), (5, 20))

    return EventTree(cn_event_tree=tree, node_to_snvs={(0, 10): {1, 12}, (5, 20): {0, 11}})


def test_ancestry_recall(event_tree: EventTree, event_tree2: EventTree):
    recall = AncestryRecall().get_scores(event_tree, event_tree2)
    assert recall['AncestryRecall'] == 1 / 4

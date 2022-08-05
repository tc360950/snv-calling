from typing import Dict, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from core.event_tree import EventTree
from core.model_data import EventTreeWithCounts
from core.types import Bin, CNEvent, EventTreeRoot


class CONETLoader:
    COUNTS_FILENAME = "inferred_counts"
    TREE_FILENAME = "inferred_tree"
    ATTACHMENT_FILENAME = "inferred_attachment"

    def __init__(self, dir: str):
        self.cc = pd.read_csv(f"{dir}cc")
        self.inferred_counts = np.loadtxt(f"{dir}{CONETLoader.COUNTS_FILENAME}").astype(
            int
        )
        self.attachment = self.__load_attachment(dir)
        self.tree = self.__load_tree(dir)

    def get_event_tree(self) -> EventTreeWithCounts:
        return EventTreeWithCounts(
            tree=EventTree(cn_event_tree=self.tree, node_to_snvs={}),
            node_to_cn_profile=dict(
                [
                    (self.attachment[cell], self.inferred_counts[cell, :])
                    for cell in range(0, len(self.attachment))
                ]
            ),
        )

    def __breakpoint_candidate_idx_to_bin(self, idx: int) -> Bin:
        counter = 0
        for i in range(0, self.cc.shape[0]):
            if self.cc.candidate_brkp.iloc[i] == 1.0:
                if counter == idx:
                    return i
                counter += 1

    def __convert_conet_node(self, node: Tuple[int, int]) -> CNEvent:
        if node != EventTreeRoot:
            return self.__breakpoint_candidate_idx_to_bin(
                node[0]
            ), self.__breakpoint_candidate_idx_to_bin(node[1])
        return EventTreeRoot

    def __load_attachment(self, dir: str) -> Dict[int, CNEvent]:
        attachment = []
        with open(f"{dir}{CONETLoader.ATTACHMENT_FILENAME}", "r") as f:
            for line in f:
                node = (int(line.split(";")[1]), int(line.split(";")[2]))
                attachment.append(self.__convert_conet_node(node))
        return dict([(i, attachment[i]) for i in range(0, len(attachment))])

    def __load_tree(self, dir: str) -> nx.DiGraph:
        tree = nx.DiGraph()

        def load_node(node: str) -> CNEvent:
            return self.__convert_conet_node(
                (int(node.split(",")[0]), int(node.split(",")[1]))
            )

        with open(f"{dir}{CONETLoader.TREE_FILENAME}", "r") as f:
            for line in f:
                line = line.replace("(", "").replace(")", "")
                tree.add_edge(
                    load_node(line.split("-")[0]), load_node(line.split("-")[1])
                )
        return tree

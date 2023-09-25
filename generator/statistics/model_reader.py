import pickle
from collections import defaultdict
from pathlib import Path

import numpy
import networkx as nx


class ModelReader:
    def __init__(self, dir):
        self.tree_path = Path(dir) / Path("event_tree")
        self.attachment_path = Path(dir) / Path("attachment")
        self.cn_path = Path(dir) / Path("cn")
        self.model_path = Path(dir) / Path("model")
        self._tree = None
        self._attachment = None
        self._cn = None
        self._snvs = None
        self._genotypes = None

    @property
    def cell_snv_pairs(self):
        result = []
        for cell, node in enumerate(self._attachment):
            path = nx.shortest_path(self.tree, (0, 0), node)
            for n in path:
                for snv in self._snvs.get(n, set()):
                    result.append((cell, snv))
        return result

    @property
    def tree(self):
        return self._tree

    @property
    def attachment(self):
        return self._attachment

    @property
    def cn(self):
        return self._cn

    @property
    def snvs(self):
        return self._snvs

    @property
    def genotypes(self):
        return self._genotypes

    @property
    def ancestor_descendant_pairs(self) -> list[tuple[int, int]]:
        result = []
        for node in self._tree.nodes:
            desc: set = nx.descendants(self._tree, node)
            if node in desc:
                desc.remove(node)
            node_cells = [c for c, n in enumerate(self._attachment) if n == node]
            desc_cells = [c for c, n in enumerate(self._attachment) if n in desc]
            for c1 in node_cells:
                for c2 in desc_cells:
                    result.append((c1, c2))
        return result

    @property
    def branching_pairs(self) -> set[tuple[int, int]]:
        result = set()
        for c1 in range(0, len(self._attachment)):
            for c2 in range(0, len(self._attachment)):
                if c1 < c2:
                    result.add((c1, c2))
        for cell_pair in self.ancestor_descendant_pairs:
            c1, c2 = min(cell_pair), max(cell_pair)
            result.remove((c1, c2))
        for node in self._tree.nodes:
            node_cells = [c for c, n in enumerate(self._attachment) if n == node]
            for c1 in node_cells:
                for c2 in node_cells:
                    if c1 < c2:
                        result.remove((c1, c2))
        return result
    def load(self):
        tree = self._load_tree()
        self._tree = tree.tree.cn_event_tree
        self._snvs = tree.tree.node_to_snvs
        self._attachment = self._load_attachment()
        self._cn = self._load_cn()
        self._genotypes = self._load_genotypes()

    def _load_genotypes(self):
        with self.model_path.open(mode="rb") as f:
            model = pickle.load(f)
        node_to_all_snvs = defaultdict(set)
        for node in model.tree.cn_event_tree.nodes:
            node_to_all_snvs[node] = model.tree.gather_snvs_on_path_to_root(node)
        genotypes = set()
        for cell, node in enumerate(self.attachment):
            for snv in node_to_all_snvs[node]:
                genotypes.add(
                    (cell, snv, int(model.node_to_altered_counts_profile[node][snv]), int(model.node_to_cn_profile[node][snv]))
                )
        return genotypes
    def _load_cn(self):
        return numpy.transpose(numpy.loadtxt(str(self.cn_path), delimiter=' ', dtype=int))

    def _load_attachment(self) -> list:
        with self.attachment_path.open(mode="rb") as f:
            return pickle.load(f)

    def _load_tree(self):
        with self.tree_path.open(mode="rb") as f:
            tree = pickle.load(f)
            return tree


if __name__ == "__main__":
    c = ModelReader("../../tmp")
    c.load()
    print("x")

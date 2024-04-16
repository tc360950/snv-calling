import pickle
from collections import defaultdict
from pathlib import Path
from typing import Tuple, List, Set

import networkx as nx
import numpy
import pandas


class ConetReader:
    def __init__(self, dir, cc_path: Path, postfix: str):
        self.cc_path = cc_path
        self.tree_path = Path(dir) / Path(f"inferred_tree")
        self.attachment_path = Path(dir) / Path(f"inferred_attachment")
        self.snvs_path = Path(dir) / Path(f"inferred_snvs2")
        self.cn_path = Path(dir) / Path(f"inferred_counts")
        self.genotypes_path = Path(dir) / Path("inferred_genotypes")
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
    def ancestor_descendant_pairs(self) -> List[Tuple[int, int]]:
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
    def not_ancestor_descendant_snv_pairs(self) -> List[Tuple[int, int]]:
        all_snvs = set()
        for _, snvs in self._snvs.items():
            all_snvs.update(snvs)
        desc = set(self.ancestor_descendant_snv_pairs)
        result = []
        for snv1 in all_snvs:
            for snv2 in all_snvs:
                if snv1 < snv2 and (snv1, snv2) not in desc and (snv2, snv1) not in desc:
                    result.append((snv1, snv2))
        return result
    @property
    def ancestor_descendant_snv_pairs(self) -> List[Tuple[int, int]]:
        result = []
        for node in self._tree.nodes:
            desc: set = nx.descendants(self._tree, node)
            if node in desc:
                desc.remove(node)
            node_snvs = [snvs for edge, snvs in self._snvs.items() if edge == node]
            if node_snvs:
                node_snvs = node_snvs[0]

            desc_snvs = set()
            for descendant in desc:
                for edge, snvs in self._snvs.items():
                    if descendant == edge:
                        desc_snvs.update(snvs)
            for n_snv in node_snvs:
                for desc_snv in desc_snvs:
                    result.append((n_snv, desc_snv))
            for n_snv1 in node_snvs:
                for n_snv2 in node_snvs:
                    if n_snv1 < n_snv2:
                        result.append((n_snv1, n_snv2))
        return result

    @property
    def branching_pairs(self) -> Set[Tuple[int, int]]:
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
        self._tree = self._load_tree()
        self._attachment = self._load_attachment()
        self._snvs = self._load_snvs()
        self._cn = self._load_cn()
        self._genotypes = self._load_genotypes()

    def _load_genotypes(self):
        x = numpy.loadtxt(str(self.genotypes_path), delimiter=';', dtype=int)
        genotypes = set()
        for i in range(0, x.shape[0]):
            genotypes.add(
                (x[i, 0], x[i,1], x[i,2], x[i,3])
            )
        return genotypes

    def _load_cn(self):
        return numpy.loadtxt(str(self.cn_path), delimiter=';', dtype=int)

    def _load_snvs(self):
        result = defaultdict(set)
        def line_to_edge(line):
            line = line.replace('1_', '')
            parent, child, rest = line.split(';')
            return int(parent), int(child), int(rest)

        with self.snvs_path.open(mode="r") as f:
            data = f.read().split('\n')
            for d in data:
                if d == '':
                    continue
                parent, child, snv = line_to_edge(d)
                result[(parent, child)].add(snv)
        return result

    def _load_attachment(self):
        with self.attachment_path.open(mode="r") as f:
            data = f.read().split('\n')
            result = []
            for d in data:
                if d == '':
                    continue
                d = d.replace('1_', '')
                cell, num, loci1, loci2 = d.split(';')
                loci1 = int(loci1)
                loci2 = int(loci2)
                if loci2 == loci1:
                    loci2 = 0
                    loci1 = 0
                result.append((loci1, loci2))
        return result

    def _load_tree(self) -> nx.DiGraph:
        with self.tree_path.open(mode="r") as f:
            data = f.read()

            def line_to_edge(line):
                line = line.replace('1_', '')
                parent, child = line.split('-')
                parent = eval(parent)
                child = eval(child)
                return parent, child

            lines = data.split('\n')
            edges = []
            for line in lines:
                if line == '':
                    continue
                edges.append(line_to_edge(line))
            graph = nx.DiGraph()
            graph.add_edges_from(edges)
            return graph


if __name__ == "__main__":
    c = ConetReader("../../tmp/out")
    c.load()
    print("x")

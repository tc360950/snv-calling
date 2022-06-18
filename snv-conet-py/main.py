import json

from generator.context import SNVGeneratorContext
from generator.model import SNVModel, CellData
from generator.gen_utils import sample_conditionally_with_replacement

from typing import List, Tuple, Optional
import numpy as np
import random

from stats.statistics_calculator import StatisticsCalculator


class GeneratorContext(SNVGeneratorContext):
    def get_max_cn(self) -> int:
        return 9

    def get_b_sampling_error(self) -> float:
        return 0.00001

    def get_d_sampling_error(self) -> float:
        return 0.01

    def get_read_success_probablity(self) -> float:
        return 0.0001

    def get_per_allele_coverage(self) -> float:
        return 0.03

    def get_neutral_cn(self) -> int:
        return 2

    def sample_cn_change(self, event: Tuple[int, int], parent_cn_profile: np.ndarray,
                         overlap_bit_map: np.ndarray) -> int:
        """
            Wyznaczanie zmiany CN przez wydarzenie @event.
            parent_cn_profile - tablica (dlugosc to ilosc binow) trzymajaca CN profile rodzica w drzewie
            overlap_bit_map - tablica (dlugosc to ilosc binow) 0, 1, gdzie oberlap_bit_map[bin] == 1.0
            oznacza, ze bin jest objety jakims wydarzeniem w potomku wierzcholka @event - a wiec nie moze
            byc na tym binie calkowitej delecji.

        """
        max_amplification = min(2, self.get_max_cn(), max([parent_cn_profile[b] for b in range(event[0], event[1])]))
        if np.sum(overlap_bit_map) == 0.0:  # No bin is present in child nodes- we can do full deletion
            return sample_conditionally_with_replacement(1, lambda: random.randint(-int(np.min(parent_cn_profile)),
                                                                                   max_amplification + 1),
                                                         lambda x: x != 0)[0]

        min_cn_in_bins_with_overlap = np.min(parent_cn_profile[overlap_bit_map == 1.0])
        max_possible_deletion = -min(2, int(min_cn_in_bins_with_overlap) - 1)
        return \
        sample_conditionally_with_replacement(1, lambda: random.randint(max_possible_deletion, max_amplification + 1),
                                              lambda x: x != 0)[0]

    def get_number_of_bins(self) -> int:
        return 100

    def get_cn_event_candidates(self) -> List[Tuple[int, int]]:
        return [(a, b) for a in range(0, self.get_number_of_bins()) for b in range(0, self.get_number_of_bins()) if a < b and b - a < 50]

    def get_snv_event_candidates(self) -> List[int]:
        """
        Lista wszystkich mozliwych SNV, mamy jeden SNV per bin, wiec SNV reprezentujemy przez nr bin
        """
        return [i for i in range(0, self.get_number_of_bins())]

    def sample_number_of_snvs_for_edge(self, num_available_snvs: int) -> int:
        """
            Losowanie ilosci SVN na krawedzi.
            @num_available_snvs = Ile jest dostepnych SNV
        """
        return min(num_available_snvs, np.random.poisson(lam=1.0, size=None))

    def get_number_of_alterations(self, cn_before: int, cn_after: int, parent_altered_counts: Optional[int]) -> int:
        if parent_altered_counts is None:  # SNV appears on the edge to the node
            if cn_after == 0:
                return 0
            if cn_after == cn_before:
                return 1
            elif cn_after > cn_before:
                return 1 if random.randint(0, cn_before) > 0 else cn_after - cn_before
            else:
                return 1 if random.randint(0, cn_before) < cn_after else 0
        else:  # SVN appeared somewhere earlier
            if cn_after == cn_before:
                return parent_altered_counts
            if cn_after > cn_before:
                return parent_altered_counts + (cn_after - cn_before) if random.randint(0,
                                                                                        cn_before) < parent_altered_counts else parent_altered_counts
            if cn_after == 0:
                return 0
            else:
                copies_for_deletion = random.sample(range(0, cn_before), cn_before - cn_after)
                return parent_altered_counts - len([x for x in copies_for_deletion if x < parent_altered_counts])


model = SNVModel(GeneratorContext())
model.generate_structure(10)



CLUSTER_SIZE = 100
CLUSTERS = 1000

from solver import CellsData, MLSolver


cluster_to_cell_data = {}

for c in range(0, CLUSTERS):
    print(f"Generating cluster num {c}")
    node = next(model.generate_cell_attachment())
    cell = model.generate_cell_in_node(node)
    for _ in range(0, CLUSTER_SIZE - 1):
        cell2 = model.generate_cell_in_node(node)
        cell.d = cell.d + cell2.d
        cell.b = cell.b + cell2.b
    cluster_to_cell_data[c] = cell


cluster_to_size = dict([(c, CLUSTER_SIZE) for c in range(0, CLUSTERS)])
attachment = dict([(c, cluster_to_cell_data[c].attachment) for c in range(0, CLUSTERS)])

b = np.vstack([cluster_to_cell_data[c].b for c in range(0, CLUSTERS)])
d = np.vstack([cluster_to_cell_data[c].d for c in range(0, CLUSTERS)])

cells_data = CellsData(d=d, b=b, attachment=attachment, cell_cluster_sizes=cluster_to_size)


tree_with_snvs = model.event_tree
model.event_tree = model.event_tree.create_copy_without_snvs()


solver = MLSolver(cells_data, model, GeneratorContext())


solver.insert_snv_events()

print("x")

stats = StatisticsCalculator.calculate_stats(model.event_tree, tree_with_snvs)

print(json.dumps(stats))
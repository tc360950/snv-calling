import pytest

import generator.gen
from core.model_data import CellsData, EventTreeWithCounts, Parameters
from core.types import EventTreeRoot
from generator.context import SNVGeneratorContext
from generator.model import SNVModel
from solver.ml_snv_solver import MLSolver


@pytest.fixture
def ctxt() -> SNVGeneratorContext:
    return SNVGeneratorContext(p=Parameters(m=0.3, q=0.0001, e=0.001))


def test_empty_model(ctxt: SNVGeneratorContext):  # Tree is just a root node - no SNV can be added
    empty_model = SNVModel.create_empty_model(ctxt)
    cells: CellsData = generator.gen.sample_cells_data(10, 100, empty_model, ctxt=ctxt, random_attachment=True)[0]
    inferred_snvs = MLSolver(
        cd=cells,
        tree=EventTreeWithCounts(tree=empty_model.tree, node_to_cn_profile=empty_model.node_to_cn_profile),
        p=Parameters(e=ctxt.sequencing_error(), m=ctxt.per_allele_coverage(), q=ctxt.read_success_prob())
    ).insert_snv_events([i for i in range(0, ctxt.number_of_bins())])
    assert inferred_snvs.node_to_snvs == {}


def test_model_without_snv(ctxt: SNVGeneratorContext):  # Tree has no SNVs
    empty_model = SNVModel.create_empty_model(ctxt)
    empty_model.tree.cn_event_tree.add_edge(EventTreeRoot, (0, 10))
    empty_model.tree.cn_event_tree.add_edge((0, 10), (5, 15))
    empty_model.node_to_cn_profile[(0, 10)] = empty_model.node_to_cn_profile[EventTreeRoot]
    empty_model.node_to_cn_profile[(5, 15)] = empty_model.node_to_cn_profile[EventTreeRoot]
    empty_model.node_to_altered_counts_profile[(0, 10)] = empty_model.node_to_altered_counts_profile[EventTreeRoot]
    empty_model.node_to_altered_counts_profile[(5, 15)] = empty_model.node_to_altered_counts_profile[EventTreeRoot]

    cells: CellsData = generator.gen.sample_cells_data(50, 100, empty_model, ctxt=ctxt, random_attachment=False)[0]
    inferred_snvs = MLSolver(
        cd=cells,
        tree=EventTreeWithCounts(tree=empty_model.tree, node_to_cn_profile=empty_model.node_to_cn_profile),
        p=Parameters(e=ctxt.sequencing_error(), m=ctxt.per_allele_coverage(), q=ctxt.read_success_prob())
    ).insert_snv_events([i for i in range(0, ctxt.number_of_bins())])

    for v in inferred_snvs.node_to_snvs.values():
        assert len(v) == 0


def test_model_with_snvs(ctxt: SNVGeneratorContext): # Tree has SNVs, sample size is sufficient - everything should be inferred correctly
    model = generator.gen.SNVModelGenerator(ctxt=ctxt).generate_model(4)
    cells: CellsData = generator.gen.sample_cells_data(50, 100, model, ctxt=ctxt, random_attachment=False)[0]
    inferred_snvs = MLSolver(
        cd=cells,
        tree=EventTreeWithCounts(tree=model.tree, node_to_cn_profile=model.node_to_cn_profile),
        p=Parameters(e=ctxt.sequencing_error(), m=ctxt.per_allele_coverage(), q=ctxt.read_success_prob())
    ).insert_snv_events([i for i in range(0, ctxt.number_of_bins())])

    assert inferred_snvs.node_to_snvs == model.tree.node_to_snvs

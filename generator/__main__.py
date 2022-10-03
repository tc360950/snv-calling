import logging
import pickle

import numpy
import numpy as np

from generator.conet_integration.utils import raw_cc_matrix_to_CONET_format
from generator.generator.context import SNVGeneratorContext
from generator.generator.gen import SNVModelGenerator, sample_cells_data

import argparse

parser = argparse.ArgumentParser(description='Generate CONET SNV tree.')
parser.add_argument('--sequencing_error', type=float,default=0.001)
parser.add_argument('--read_success_prob', type=float, default=0.0001)
parser.add_argument('--coverage', type=float,default=0.03)
parser.add_argument('--max_cn', type=int, default=9)
parser.add_argument('--neutral_cn', type=int, default=2)
parser.add_argument('--bins', type=int, default=1000)
parser.add_argument('--tree_size', type=int, default=20)
parser.add_argument('--clusters', type=int, default=20, help="Number of cell clusters")
parser.add_argument('--cluster_size', type=int, default=20)
parser.add_argument('--random_attachment', type=bool, default=True, help="If set to False clusters will be attached to node in deterministic fashion")
parser.add_argument('--snv_event_proportion', type=float, default=0.5,
                    help="average proportion of SNV events on the tree")

args = parser.parse_args()


class Context(SNVGeneratorContext):

    def __init__(self):
        pass

    def sequencing_error(self) -> float:
        return args.sequencing_error

    def read_success_prob(self) -> float:
        return args.read_success_prob

    def per_allele_coverage(self) -> float:
        return args.coverage

    def max_cn(self) -> int:
        return args.max_cn

    def neutral_cn(self) -> int:
        return args.neutral_cn

    def number_of_bins(self) -> int:
        return args.bins

    def snv_event_prob(self):
        return args.snv_event_proportion


if __name__ == "__main__":
    ctxt = Context()
    model_generator = SNVModelGenerator(ctxt)
    model = model_generator.generate_model(args.tree_size)

    while not any([np.min(n) == 0.0 for n in model.node_to_cn_profile.values()]):
        logging.getLogger().error("Failed to generate tree with at least one deletion, trying one more time...")
        model = model_generator.generate_model(args.tree_size)

    cells_data, cc_matrix = sample_cells_data(clusters=args.clusters, cluster_size=args.cluster_size, model=model,
                                              ctxt=ctxt, random_attachment=args.random_attachment)

    raw_cc_matrix_to_CONET_format(cc_matrix, model.tree).to_csv(f"./cc", index=False)
    with open("./cluster_sizes", "w") as f:
        for _ in range(0, args.clusters):
            f.write(f"{args.cluster_size}\n")

    with open("./snvs_data", "w") as f:
        breakpoints = model.tree.get_breakpoint_loci()
        for i in range(0, ctxt.number_of_bins()):
            smaller = [b for b in breakpoints if b <= i]
            if smaller:
                f.write(f"{i};{len(smaller) - 1}\n")
            else:
                f.write(f"{i};{-1}\n")

    numpy.savetxt("B", cells_data.b, delimiter=";")
    numpy.savetxt("D", cells_data.d, delimiter=";")
    with open(f"./event_tree", "wb") as f:
        pickle.dump(model.tree.event_tree, f)

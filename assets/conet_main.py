import argparse

import numpy as np
import pandas as pd

from conet.data_converter.corrected_counts import CorrectedCounts
from conet.data_converter.data_converter import DataConverter
from conet import CONET, CONETParameters, InferenceResult

parser = argparse.ArgumentParser(description='Run CONET')
parser.add_argument('--corrected_counts_file', type=str, required=True)
parser.add_argument('--data_dir', type=str, default='./')
parser.add_argument('--param_inf_iters', type=int, default=100000)
parser.add_argument('--pt_inf_iters', type=int, default=100000)
parser.add_argument('--counts_penalty_s1', type=float, default=0.0)
parser.add_argument('--counts_penalty_s2', type=float, default=0.0)
parser.add_argument('--event_length_penalty_k0', type=float, default=1.0)
parser.add_argument('--tree_structure_prior_k1', type=float, default=1.0)
parser.add_argument('--use_event_lengths_in_attachment', type=bool, default=True)
parser.add_argument('--seed', type=int, default=12312)
parser.add_argument('--mixture_size', type=int, default=4)
parser.add_argument('--num_replicas', type=int, default=5)
parser.add_argument('--threads_likelihood', type=int, default=4)
parser.add_argument('--verbose', type=bool, default=True)
parser.add_argument('--neutral_cn', type=float, default=2.0)
parser.add_argument('--output_dir', type=str, default='./')

args = parser.parse_args()

if __name__ == "__main__":
    corrected_counts: pd.DataFrame = pd.read_csv(args.corrected_counts_file)
    cc = CorrectedCounts(corrected_counts)
    DataConverter(event_length_normalizer=corrected_counts.shape[0]).create_CoNET_input_files(out_path=args.data_dir,
                                                                                              corrected_counts=cc)
    conet = CONET("./CONET")
    params = CONETParameters(
        data_dir=args.data_dir,
        param_inf_iters=args.param_inf_iters,
        pt_inf_iters=args.pt_inf_iters,
        counts_penalty_s1=args.counts_penalty_s1,
        counts_penalty_s2=args.counts_penalty_s2,
        event_length_penalty_k0=args.event_length_penalty_k0,
        tree_structure_prior_k1=args.tree_structure_prior_k1,
        use_event_lengths_in_attachment=args.use_event_lengths_in_attachment,
        seed=args.seed,
        mixture_size=args.mixture_size,
        num_replicas=args.num_replicas,
        threads_likelihood=args.threads_likelihood,
        verbose=args.verbose,
        neutral_cn=args.neutral_cn,
        output_dir=args.output_dir
    )
    conet.infer_tree(params)
    result = InferenceResult(params.output_dir, cc)
    inferred_counts = np.transpose(result.get_inferred_copy_numbers(neutral_cn=int(params.neutral_cn)))
    np.savetxt(f"{params.output_dir}inferred_counts", inferred_counts)

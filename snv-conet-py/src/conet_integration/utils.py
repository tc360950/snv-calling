import numpy as np
import pandas as pd

from core.event_tree import EventTree


def raw_cc_matrix_to_CONET_format(cc: np.ndarray, event_tree: EventTree) -> pd.DataFrame:
    no_cells, no_loci = cc.shape
    cc = np.transpose(cc)

    # Add columns required by CONET
    add = np.full([no_loci, 5], 1.0, dtype=np.float64)
    add[:, 1] = range(0, no_loci)  # Bin start
    add[:, 2] = range(1, no_loci + 1)  # Bin end
    add[:, 4] = 0  # Breakpoint markers
    add[event_tree.get_breakpoint_loci(), 4] = 1
    full_counts = np.hstack([add, cc])
    return pd.DataFrame(full_counts,
                        columns=["chr", "start", "end", "width", "candidate_brkp"] + ["cell" + str(i) for i in
                                                                                      range(0, no_cells)])
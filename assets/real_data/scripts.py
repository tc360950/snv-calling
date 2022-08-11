import pandas as pd
from conet.data_converter.corrected_counts import CorrectedCounts

data = pd.read_csv("./cn_matrix.csv", sep=",", header=0)
CC = pd.read_csv("./cc_matrix.csv", sep=",", header=0)
CC["candidate_brkp"] = 0
CC["width"] = 150000

candidates = []
for i in range(1, data.shape[0]):
    if len(set(data.iloc[i, 4:]).difference(set(data.iloc[i - 1, 4:]))) > 0:
        if not any(CC.iloc[i].isnull()):
            candidates.append(i)

CC.candidate_brkp.iloc[candidates] = 1

CC.chr = CC.chr.apply(lambda x: 23 if x == "X" else 24 if x == "Y" else int(x))
CC.fillna(2.0, inplace=True)

cols = CC.columns.tolist()
cols = cols[1 : (len(cols) - 2)]
cols.insert(3, "width")
cols.insert(4, "candidate_brkp")
CC = CC[cols]

# For this install conet-py with pip install
# git+https://github.com/tc360950/CONET/#subdirectory=python/conet-py
cc_ = CorrectedCounts(CC)
cc_.add_chromosome_ends(neutral_cn=2, end_bin_length=150000)
CC = cc_.corr_counts_df
CC.to_csv("cc", sep=",", index=False)

snv_to_bin = pd.read_csv("snv_to_bin_dict.csv", sep=",", header=0)
snv_to_bin["chr"] = snv_to_bin.snvs.apply(lambda x: x.split("_")[0])
snv_to_bin["pos"] = snv_to_bin.snvs.apply(lambda x: int(x.split("_")[1]))
snv_to_bin.chr = snv_to_bin.chr.apply(
    lambda x: 23 if x == "X" else 24 if x == "Y" else int(x)
)
snv_to_bin.sort_values(by=["chr", "pos"], inplace=True)
for i in range(0, snv_to_bin.shape[0]):
    snv_to_bin.python_row.iloc[i] += snv_to_bin.chr.iloc[i] - 1

snv_to_bin = snv_to_bin[["python_row"]]
snv_to_bin.to_csv("snv_to_bin", header=False, index=False)

d = pd.read_csv("d_matrix.csv", header=0)
b = pd.read_csv("b_matrix.csv", header=0)

d["chr"] = d.snvs.apply(lambda x: x.split("_")[0])
d.chr = d.chr.apply(lambda x: 23 if x == "X" else 24 if x == "Y" else int(x))
d["pos"] = d.snvs.apply(lambda x: int(x.split("_")[1]))

b["chr"] = b.snvs.apply(lambda x: x.split("_")[0])
b.chr = b.chr.apply(lambda x: 23 if x == "X" else 24 if x == "Y" else int(x))
b["pos"] = b.snvs.apply(lambda x: int(x.split("_")[1]))

d.sort_values(by=["chr", "pos"], inplace=True)
b.sort_values(by=["chr", "pos"], inplace=True)

b = b.drop(["chr", "pos", "snvs", "Unnamed: 0"], axis=1)
d = d.drop(["chr", "pos", "snvs", "Unnamed: 0"], axis=1)
d.to_csv("d", index=False, sep=";")
b.to_csv("b", index=False, sep=";")
if __name__ == "__main__":
    print("x")

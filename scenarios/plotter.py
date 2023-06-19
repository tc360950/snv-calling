import matplotlib.pyplot as plt
import pandas as pd

file_path = './scenarios/scenario2/results_1'
df = pd.read_csv(file_path, delimiter=";")

print(df.head())

print(df.median())

df = df[df.scenario.isin(["tree_unknown_brkp_known_0.1", "tree_known_brkp_known", "tree_known_attachment_unknown_brkp_known", "tree_unknown_brkp_unknown_0.1"])]

fig, ax = plt.subplots()
medians = df.groupby(['scenario'])['cell_snv_recall'].median()

boxplot = df.boxplot(column='cell_snv_precision', by='scenario', ax=ax)
positions = range(1, len(medians) + 1)
for tick, value in zip(positions, medians):
    boxplot.axhline(value, xmin=(tick-1)/len(medians), xmax=tick/len(medians), color='r', linestyle='dashed', linewidth=2)
plt.xticks(rotation=25)
plt.show()
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file into a DataFrame
file_path = '../wyniki'  # Replace with your file path
df = pd.read_csv(file_path, delimiter=";")

# Check the DataFrame
print(df.head())
#scenario;inferred_tree_size;real_tree_size;inferred_snvs;real_snvs;cn_node_recall;cn_node_precision;cn_edge_recall;cn_edge_precision;snv_recall;snv_precision;cell_snv_recall;cell_snv_precision;cn_prob;non_trivial_cns

MEASURE = 'snv_recall'
# Create a box plot of values per each scenario value
fig, ax = plt.subplots()
df.boxplot(column=MEASURE, by='scenario', ax=ax)

# Customize the plot
ax.set_title(f'Box Plot of {MEASURE} per Scenario')
ax.set_xlabel('Scenario')
ax.set_ylabel(MEASURE)

# Show the plot
plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file into a DataFrame
file_path = './results_test3'  # Replace with your file path
df = pd.read_csv(file_path, delimiter=";")

# Check the DataFrame
print(df.head())

print(df.median())

import pandas as pd
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
# Assuming df is your DataFrame
medians = df.groupby(['scenario'])['cell_snv_recall'].median()

# Create boxplot
boxplot = df.boxplot(column='cell_snv_precision', by='scenario', ax=ax)

# Get the x position of each box
positions = range(1, len(medians) + 1)

# Add median lines over each box
for tick, value in zip(positions, medians):
    boxplot.axhline(value, xmin=(tick-1)/len(medians), xmax=tick/len(medians), color='r', linestyle='dashed', linewidth=2)
plt.xticks(rotation=25)
plt.show()
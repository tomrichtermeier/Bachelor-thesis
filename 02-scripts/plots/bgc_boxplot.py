import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/data/FUNC_antiSMASH.tsv", sep='\t')
df = pd.DataFrame(data)

def group_by_sample_id(sample_id):
    if 'UT' in sample_id:
        return 'UT'
    elif 'AW' in sample_id:
        return 'AW'
    elif 'Zape' in sample_id:
        return 'Zape'
    
# Group individual samples to region (AW, UT, Zape) and calculate abundances
df['Sample_Group'] = df['Sample_ID'].apply(group_by_sample_id)
result = df.groupby(['Sample_Group', 'Product_class']).size().reset_index(name='Abundance')
result= pd.DataFrame(result)

boxplot_data = df.iloc[:,[8,14]]



df = pd.DataFrame(data)

fig, ax = plt.subplots(figsize=(6, 6))

order = ['Zape', 'AW', 'UT']
colors = {'Zape': 'cornflowerblue', 'UT': 'lightgreen', 'AW': 'lightcoral'}

boxplot_data = [df[df['Sample_Group'] == group]['BGC_length'] for group in order]

# Create the boxplot
boxplot = ax.boxplot(boxplot_data, labels=order)

# Set colors for data points (outliers)
for i, group in enumerate(order):
    outliers = boxplot['fliers'][i]
    outliers.set(marker='o', markersize=6, markeredgecolor='black', markerfacecolor=colors[group], alpha=0.6)
for i, group in enumerate(order):
    median = boxplot['medians'][i]
    median.set_color(colors[group])
    
ax.set_ylabel('BGC length [bp]')
ax.grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax.spines[['top', 'right']].set_visible(False)

plt.show()

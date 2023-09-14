import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_results/data/FUNC_antiSMASH.tsv", sep='\t')
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


fig, ax = plt.subplots(1,2, figsize=(6.5, 4), gridspec_kw={
                       'width_ratios': [1.7, 1]})

order = ['Zape', 'AW', 'UT']
colors = {'Zape': 'cornflowerblue', 'UT': 'lightgreen', 'AW': 'lightcoral'}

boxplot_data = [df[df['Sample_Group'] == group]['BGC_length'] for group in order]

# Create the boxplot
boxplot = ax[0].boxplot(boxplot_data, labels=order)

# Set colors for data points (outliers)
for i, group in enumerate(order):
    outliers = boxplot['fliers'][i]
    outliers.set(marker='o', markersize=6, markeredgecolor='black', markerfacecolor=colors[group], alpha=0.6)
#for i, group in enumerate(order):
 #   median = boxplot['medians'][i]
  #  median.set_color(colors[group])
ax[0].set_title('a', loc='left', fontweight="bold")
ax[0].set_ylabel('BGC length [bp]')
ax[0].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax[0].spines[['top', 'right']].set_visible(False)




data = pd.read_csv(
    "/Users/tomrichtermeier/Documents/bachelor_results/data/FUNC_antiSMASH.tsv", sep='\t')
df = pd.DataFrame(data)

complete_contigs = df.iloc[:, [0, 2, 3, 5, 8]]
complete_contigs = complete_contigs[complete_contigs['BGC_complete'] == 'Yes']
# create boxplot for length distribution
flierprops = dict(marker='o', markerfacecolor='silver', markersize=6, linestyle='none', markeredgecolor='gray', alpha=0.6)
#median = median.set_color('lightcoral'
ax[1].boxplot(complete_contigs.iloc[:, 4], flierprops=flierprops, medianprops={'color': 'orange'})
ax[1].set_ylabel('BGC length [bp]')
ax[1].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax[1].spines[['top', 'right']].set_visible(False)
ax[1].set_xticks([1])
ax[1].set_xticklabels([''])
ax[1].set_title('b', loc='left', fontweight="bold")
#ax[1].boxplot(complete_contigs.iloc[:, 4], flierprops=flierprops)
plt.tight_layout()
plt.show()
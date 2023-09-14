import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv(
    "/Users/tomrichtermeier/Documents/bachelor_results/data/FUNC_antiSMASH.tsv", sep='\t')
df = pd.DataFrame(data)

complete_contigs = df.iloc[:, [0, 2, 3, 5, 8]]
complete_contigs = complete_contigs[complete_contigs['BGC_complete'] == 'Yes']

bgc_types = complete_contigs['Product_class'].value_counts()

# PLOTS
fig, ax = plt.subplots(figsize=(6, 4)) #figsize=(6, 6)

# create barplot for bgc types
categories = bgc_types.index[0:9]
categories_short = ['RRE-containing', 'Ranthipeptide', 'CLA', 'NRPS-like',  'Arylpolyene', 'Betalactone',
                    'RiPP-like', 'Terpene', 'NRPS']  # names are too long for figure
x = np.arange(len(categories))
ax.bar(x, bgc_types.iloc[0:9], width=0.33,
          color='cornflowerblue', zorder=3)
ax.grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax.spines[['top', 'right']].set_visible(False)
ax.set_xticks(x)
ax.set_xticklabels(categories_short, rotation=90)
ax.set_ylabel('Amount of complete BGCs across all sampels')

plt.tight_layout()
plt.show()
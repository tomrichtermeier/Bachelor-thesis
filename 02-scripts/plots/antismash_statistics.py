#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 14:12:27 2023

@author: tomrichtermeier
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
data = pd.read_csv(
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/FUNC_antiSMASH.tsv", sep='\t')
df = pd.DataFrame(data)

# General statistics
number_regions = df.iloc[:, 0].value_counts()
completeness = df.iloc[:, [0, 5]].value_counts()
avg_size = df.iloc[:, [0, 8]].value_counts()
productclasses_total = df['Product_class'].value_counts()

completeness = df.iloc[:, [0, 8]]
completeness_avg = completeness.groupby(
    'Sample_ID')['BGC_length'].mean().reset_index().round()

# BGC-Classes of all samples combined
bgc_types = df['Product_class'].value_counts()
bgc_types = pd.DataFrame(bgc_types)


fig, ax = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={
                       'width_ratios': [1, 1.6]})


# BGC-Classes of AW,UT and Zape samples
def group_by_sample_id(sample_id):
    if 'UT' in sample_id:
        return 'UT'
    elif 'AW' in sample_id:
        return 'AW'
    elif 'Zape' in sample_id:
        return 'Zape'


# Group individual samples to region (AW, UT, Zape) and calculate abundances
df['Sample_Group'] = df['Sample_ID'].apply(group_by_sample_id)
result = df.groupby(['Sample_Group', 'Product_class']
                    ).size().reset_index(name='Abundance')
result = pd.DataFrame(result)

# sort
result = result.sort_values(
    ['Sample_Group', 'Abundance'], ascending=[True, False])

# Transform abundances to percantage
total_number = result.groupby('Sample_Group')['Abundance'].transform('sum')
result['Relative_Abundance'] = (result['Abundance']/total_number)*100


# FIGURE
def get_relabundance(Sample_Group, Product_class):
    abundance = result.loc[(result['Sample_Group'] == Sample_Group) & (
        result['Product_class'] == Product_class), 'Relative_Abundance']
    return abundance.values[0]


categories = ['Cyclic-lactone-autoinducer', 'RRE-containing', 'Ranthipeptide',
              'Terpene', 'Arylpolyene', 'NRPS', 'RiPP-like', 'NRPS-like', 'Phosphonate']
zape = [get_relabundance('Zape', categories[0]),
        get_relabundance('Zape', categories[1]),
        get_relabundance('Zape', categories[2]),
        get_relabundance('Zape', categories[3]),
        get_relabundance('Zape', categories[4]),
        get_relabundance('Zape', categories[5]),
        get_relabundance('Zape', categories[6]),
        get_relabundance('Zape', categories[7]),
        get_relabundance('Zape', categories[8])]

aw = [get_relabundance('AW', categories[0]),
      get_relabundance('AW', categories[1]),
      get_relabundance('AW', categories[2]),
      get_relabundance('AW', categories[3]),
      get_relabundance('AW', categories[4]),
      get_relabundance('AW', categories[5]),
      get_relabundance('AW', categories[6]),
      get_relabundance('AW', categories[7]),
      get_relabundance('AW', categories[8])]

ut = [get_relabundance('UT', categories[0]),
      get_relabundance('UT', categories[1]),
      get_relabundance('UT', categories[2]),
      get_relabundance('UT', categories[3]),
      get_relabundance('UT', categories[4]),
      get_relabundance('UT', categories[5]),
      get_relabundance('UT', categories[6]),
      get_relabundance('UT', categories[7]),
      get_relabundance('UT', categories[8])]

# Breite der Balken
bar_width = 0.2

# Positions
x = np.arange(len(categories))
x_tick_positions = x + bar_width
categories_short = ['CLA', 'RRE-containing', 'Ranthipeptide', 'Terpene', 'Arylpolyene',
                    'NRPS', 'RiPP-like', 'NRPS-like', 'Phosphonate']  # names are too long for figure

ax[0].bar(x, bgc_types.iloc[0:9, 0], width=0.33,
          color='cornflowerblue', zorder=3)
ax[0].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax[0].spines[['top', 'right']].set_visible(False)
ax[0].set_xticks(x)
ax[0].set_xticklabels(categories_short, rotation=90)
ax[0].set_ylabel('Amount of BGCs across all sampels')
ax[0].set_title('a', loc='left', fontweight="bold")

# Balken erstellen
ax[1].bar(x - bar_width, zape, bar_width, label='Zape',
          color='cornflowerblue', zorder=3)
ax[1].bar(x, aw, bar_width, label='AW', color='lightcoral', zorder=3)
ax[1].bar(x + bar_width, ut, bar_width,
          label='UT', color='lightgreen', zorder=3)


ax[1].set_xticks(x)
ax[1].set_title('b', loc='left', fontweight="bold")
ax[1].set_xticklabels(categories_short, rotation=90)
ax[1].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax[1].spines[['top', 'right']].set_visible(False)
ax[1].set_ylabel('Percent %')

plt.legend()
plt.rcParams['figure.dpi'] = 150
plt.tight_layout()
plt.show()


# number of binned contigs
samples = ["AW107", "AW108", "AW110", "UT30",
           "UT43", "Zape1", "Zape2", "Zape3"]
file_paths = [
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/aw107_metawrap_50_10_bins.contigs",
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/aw108_metawrap_50_10_bins.contigs",
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/aw110_metawrap_50_10_bins.contigs",
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/ut30_metawrap_50_10_bins.contigs",
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/ut43_metawrap_50_10_bins.contigs",
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/zape1_metawrap_50_10_bins.contigs",
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/zape2_metawrap_50_10_bins.contigs",
    "/Users/tomrichtermeier/Documents/bachelor_scripts/data/bins/aw107_metawrap_50_10_bins.contigs"
]

binned_contigs = {'Sample': samples, 'binned': [0]*len(samples)}
binned_contigs = pd.DataFrame(binned_contigs)

for i, file_path in enumerate(file_paths):
    data = pd.read_csv(file_path, delimiter="\t",
                       header=None,  names=['contig', 'bin'])
    contigs = set(data['contig'].tolist())
    number = df['Contig_ID'].isin(contigs).sum()
    binned_contigs.at[i, 'binned'] = number

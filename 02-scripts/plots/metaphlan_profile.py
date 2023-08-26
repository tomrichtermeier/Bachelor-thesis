#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 16:11:07 2023

@author: tomrichtermeier
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv("/mnt/archgen/users/richtermeier/bachelor_thesis/05-results/COMP_MetaPhlAn4_species_profiles.tsv", sep='\t')
df = pd.DataFrame(data)

columns=['clade_name', 'relative_abundance']
genus = df[columns].copy()
genus['clade_name'] = genus['clade_name'].str.extract(r'p__(.*?)\|')
genus_avg = genus.groupby('clade_name')['relative_abundance'].mean().reset_index()
genus_avg= genus_avg.sort_values('relative_abundance', ascending=False)

genus_avg.at[2, 'clade_name'] = 'unclassified'



fig, ax = plt.subplots()
ax.barh(genus_avg.iloc[0:17,0], genus_avg.iloc[0:17,1], color='cornflowerblue', zorder=3)
ax.set_xlabel("relative abundance [%]")
ax.grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax.spines[['top','right']].set_visible(False)


plt.rcParams['figure.dpi']=500
plt.show()

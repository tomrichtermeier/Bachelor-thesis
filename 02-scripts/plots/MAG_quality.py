#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 12:26:08 2023

@author: tomrichtermeier
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


quality={'sample': ["Zape1", "Zape2", "Zape3","AW107", "AW108", "AW110", "UT30", "UT43"],
         'HQ': [0,0,0,0,0,0,0,0],
         'LQ': [0,0,0,0,0,0,0,0]}
df_quality = pd.DataFrame(quality)



data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/Zape1-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[0,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[0,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[0,1]

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/Zape2-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[1,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[1,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[1,1]

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/Zape3-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[2,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[2,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[2,1]

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/AW107-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[3,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[3,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[3,1]

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/AW108-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[4,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[4,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[4,1]

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/AW110A-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[5,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[5,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[5,1]

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/UT30-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[6,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[6,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[6,1]

data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/binning quality/UT43-megahit_metawrap.stats", delimiter="\t")
df = pd.DataFrame(data)
df_quality.iloc[7,1] = ((df['completeness'] > 90) & (df['contamination'] < 5)).sum()
df_quality.iloc[7,2] = ((df['completeness'] >= 50) & (df['contamination'] < 10)).sum()-df_quality.iloc[7,1]




fig, ax = plt.subplots(1, 2, figsize=(11, 5))
bars = ("Zape1", "Zape2", "Zape3","AW107", "AW108", "AW110", "UT30", "UT43")
x_pos = np.arange(len(bars))
bar_width = 0.8

ax[0].bar(x_pos, df_quality.iloc[:,1], width=bar_width, label='High Quality MAG', color='lightcoral', zorder=3)
ax[0].bar(x_pos, df_quality.iloc[:,2], width=bar_width, bottom=df_quality.iloc[:,1], label='Low Quality MAG', color='cornflowerblue', zorder=3)
ax[0].set_ylabel("Amount of MAGs")
ax[0].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax[0].spines[['top','right']].set_visible(False)

ax[0].set_xticks(x_pos)
ax[0].set_xticklabels(bars, rotation=45)
# ax[0].set_yticks([10, 30, 50, 100, 150, 200, 250, max(df_quality.iloc[:,1]+df_quality.iloc[:,2])])
# ax[0].set_yticklabels(['10', '30', '50', '100', '150', '200', '250', max(df_quality.iloc[:,1]+df_quality.iloc[:,2])])
ax[0].set_yticks([10, 30, 50, 100, 150, 200, 250, 300])
ax[0].set_yticklabels(['10', '30', '50', '100', '150', '200', '250', '300'])
ax[0].set_title('a',loc='left', fontweight="bold")
ax[0].legend()


data = pd.read_csv("/mnt/archgen/users/richtermeier/bachelor_thesis/05-results/gtdbtk.summary.tsv", sep='\t')
df = pd.DataFrame(data)

phylum = df['classification'].str.extract(r'g__(.*?);')
phylum = phylum.iloc[:,0].value_counts()
phylum.rename('abundance', inplace=True)
phylum = pd.DataFrame(phylum)
phylum.index.values[1] = 'unclassified'


ax[1].barh(phylum.index[0:20], phylum.iloc[0:20,0], color='cornflowerblue', zorder=3)
ax[1].set_xlabel("Amount of classified MAGs")
ax[1].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax[1].spines[['top','right']].set_visible(False)
ax[1].set_title('b',loc='left', fontweight="bold")
# ax[1].set_xticks([0, 25, 50, 100, 150, 200, 250, 300, 350, 400])
# ax[1].set_xticklabels(['0', '25', '50', '100', '150', '200', '250', '300', '350', '400'])
#ax[1].set_xscale('log')

###### 23 cases, wo phylum reihe mit leerem namen enthält -> 23 nicht auf dem level classified
###### 500 others




plt.rcParams['figure.dpi']=300
plt.tight_layout()
plt.show()
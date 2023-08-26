#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:41:10 2023

@author: tomrichtermeier
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


### READ LENGTH DISTRIBUTION
data = pd.read_csv("/Users/tomrichtermeier/Documents/bachelor_scripts/data/fastqc_sequence_length_distribution_plot-2-1.tsv", sep='\t')
reihenfolge = ['Sequence Length (bp)', 'SRR12557706_1', 'SRR12557705_1', 'SRR12557704_1', 'SRR12557722_1', 'SRR12557711_1', 'SRR12557707_1', 'SRR12557734_1', 'SRR12557733_1']
data = data.reindex(columns=reihenfolge)
df = pd.DataFrame(data, index=['Readlänge', 'Zape1', 'Zape2', 'Zape3', 'AW107', 'AW108', 'AW110', 'UT30', 'UT43'])
fig, axs = plt.subplots(1, 1, figsize=(7.5, 5))


for i in range(1,9):
    axs.plot(data.iloc[:,0], data.iloc[:,i], label=df.index[i])
axs.set_yscale("log") 
axs.grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs.set_ylabel("number of reads")
axs.set_xlabel("read length [bp]")
axs.legend()
axs.set_title('a',loc='left', fontweight="bold")
axs.spines[['right', 'top']].set_visible(False)





# data_list = [105,95,103,90,90,105,123,128]
# zape_length = [105,95,103]
# AW1XX_length = [90,90,105]
# UT_length = [123,128]
# data_box = np.array(data_list)

# axs[1].boxplot(data_box, vert=False)
# jitter = np.random.normal(loc=0, scale=0.2, size=[3])


# axs[1].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
# axs[1].spines[['top','right']].set_visible(False)
# axs[1].plot(zape_length + jitter, np.ones(3), 'ro', alpha=0.5, label='Zape')
# axs[1].plot(AW1XX_length + jitter, np.ones(3), 'bX', alpha=0.5, label='AW')
# axs[1].plot(UT_length, np.ones(2), 'gs', alpha=0.5, label='UT')
# axs[1].set_title('b',loc='left', fontweight="bold")
# axs[1].set_xlabel("average read length [bp]")
# axs[1].set_yticks([])
# axs[1].set_yticklabels([])




plt.rcParams['figure.dpi']=100
#plt.tight_layout()

plt.show()
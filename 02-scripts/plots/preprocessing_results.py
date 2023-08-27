#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 12:00:59 2023

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
fig, axs = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={'width_ratios': [1.8, 1]})


for i in range(1,9):
    axs[0].plot(data.iloc[:,0], data.iloc[:,i], label=df.index[i])
axs[0].set_yscale("log") 
axs[0].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs[0].set_ylabel("number of reads")
axs[0].set_xlabel("read length [bp]")
axs[0].legend()
axs[0].set_title('a',loc='left', fontweight="bold")
axs[0].spines[['right', 'top']].set_visible(False)

bars = ("Zape3", "Zape2", "Zape1","AW110A", "AW108", "AW107", "UT43.2", "UT30.3")
input_reads = [138248288, 141150064, 169853426, 124607641,152273202,108548593,182274727,459767358]
output_reads = [138202346, 141025581, 169786684, 124498642, 152169729, 108507947,182253824,459695896]
diff_percent = 100-(np.array(output_reads)/np.array(input_reads))*100
zape_diff = np.array(diff_percent[0:3])
aw_diff = np.array(diff_percent[3:6])
ut_diff = np.array(diff_percent[6:9])




x1 = np.ones_like(zape_diff) * 1.22
axs[1].plot(x1,zape_diff, color='cornflowerblue', label='Zape', marker='o', linestyle = 'dotted', linewidth=0.25)


x2 = np.ones_like(aw_diff) * 1.5  # x-Werte für Gruppe 2
axs[1].plot(x2, aw_diff, color='lightcoral', label='AW', marker='s', linestyle = 'dotted', linewidth=0.25)

x3 = np.ones_like(ut_diff) * 1.78  # x-Werte für Gruppe 2
axs[1].plot(x3, ut_diff, color='lightgreen', label='UT', marker='D', linestyle = 'dotted', linewidth=0.25)


bars = ("Zape3", "Zape2", "Zape1","AW110A", "AW108", "AW107", "UT43.2", "UT30.3")
y_pos = np.arange(len(bars))
#bar1=axs[1][1].plot(y_pos, avg_read, zorder=3, linewidth=0.7)
axs[1].set_ylabel("% of reads filtered out")
#axs[1].set_title('b',loc='left', fontweight="bold")
axs[1].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs[1].spines[['top','right']].set_visible(False)
#axs[1][1].bar_label(bar1, padding=3, color='black', fontsize=8)
axs[1].set_xticks([1.1, 1.22, 1.5, 1.78, 1.9])
axs[1].set_xticklabels(['','Zape', 'AW', 'UT',''])
axs[1].set_title('b',loc='left', fontweight="bold")



plt.rcParams['figure.dpi']=200
plt.show()

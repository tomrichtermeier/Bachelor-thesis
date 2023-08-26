#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 08:25:33 2023

@author: tomrichtermeier
"""

import numpy as np
import matplotlib.pyplot as plt


bars = ("Zape1", "Zape2", "Zape3","AW107", "AW108", "AW110", "UT30", "UT43")
total_length = [606726477, 368951794, 385910544, 498460808, 442371646, 367924578, 1695566044, 547759373]
largest_contig = [144417, 176282, 268498, 144459, 224180, 135078, 461380, 192297]
contigs_10 = [4749,3217,3322,2525,3256,3555,24401,5132]
contigs_50 = [107,72,240,56,251,58,2232,201]
y_pos = np.arange(len(bars))
total_length=np.array(total_length)/10**9

fig, axs = plt.subplots(2, 2, figsize=(12, 9))  
colors = ['cornflowerblue' if 'Zape' in bar else 'lightcoral' if 'AW' in bar else 'lightgreen' for bar in bars]

axs[0,0].bar(y_pos, total_length, zorder=3, linewidth=0.7, color=colors)
axs[0,0].set_ylabel("Total length [1e9 bp]")
axs[0,0].set_title('a',loc='left', fontweight="bold")
axs[0,0].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs[0,0].spines[['top','right']].set_visible(False)


axs[0,1].bar(y_pos, largest_contig, zorder=3, linewidth=0.7, color=colors)
axs[0,1].set_ylabel("Length of largest contig [kb]")
axs[0,1].set_title('b', loc='left', fontweight="bold")
axs[0,1].set_yticks([100000, 200000,300000, 400000, 500000])
axs[0,1].set_yticklabels(['100', '200', '300', '400', '500'])
axs[0,1].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs[0,1].spines[['top','right']].set_visible(False)


axs[1,0].bar(y_pos, contigs_10, zorder=3, linewidth=0.7, color=colors)
axs[1,0].set_ylabel("Number of contigs \u2265 10 kb")
axs[1,0].set_title('c',loc='left', fontweight="bold")
axs[1,0].set_yscale('log')
axs[1,0].set_yticks([1000, 2000,3000, 4000, 5000, 6000, 10000, 20000, 25000])
axs[1,0].set_yticklabels(['1000', '2000', '3000', '4000', '5000', '6000', '10000', '20000', '25000'])
axs[1,0].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs[1,0].spines[['top','right']].set_visible(False)


axs[1,1].bar(y_pos, contigs_50, zorder=3, linewidth=0.7, color=colors)
axs[1,1].set_ylabel("Number of contigs \u2265 50 kb")
axs[1,1].set_title('d', loc='left', fontweight="bold")
axs[1,1].set_yscale('log')
axs[1,1].set_yticks([10, 50, 100, 200, 300, 500, 1000, 2000])
axs[1,1].set_yticklabels(['10', '50', '100', '200', '300', '500', '1000', '2000'])
axs[1,1].grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs[1,1].spines[['top','right']].set_visible(False)


axs[0,0].set_xticks(y_pos)
axs[0,1].set_xticks(y_pos)
axs[1,0].set_xticks(y_pos)
axs[1,1].set_xticks(y_pos)
axs[0,0].set_xticklabels(bars, rotation=45)
axs[1,0].set_xticklabels(bars, rotation=45)
axs[1,1].set_xticklabels(bars, rotation=45)
axs[0,1].set_xticklabels(bars, rotation=45)

plt.tight_layout()
plt.rcParams['figure.dpi']=500
plt.show()
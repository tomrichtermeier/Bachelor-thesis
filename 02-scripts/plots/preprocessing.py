#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:21:53 2023

@author: tomrichtermeier
"""

import numpy as np
import matplotlib.pyplot as plt

bars = ("Zape3", "Zape2", "Zape1","AW110A", "AW108", "AW107", "UT43.2", "UT30.3")
input_reads = [138248288, 141150064, 169853426, 124607641,152273202,108548593,182274727,459767358]
output_reads = [138202346, 141025581, 169786684, 124498642, 152169729, 108507947,182253824,459695896]
diff_percent = 100-(np.array(output_reads)/np.array(input_reads))*100
zape_diff = np.array(diff_percent[0:3])
aw_diff = np.array(diff_percent[3:6])
ut_diff = np.array(diff_percent[6:9])



fig, axs = plt.subplots(1, 1, figsize=(4, 4))

x1 = np.ones_like(zape_diff) * 1.22
axs.plot(x1,zape_diff, color='cornflowerblue', label='Zape', marker='o', linestyle = 'dotted', linewidth=0.25)


x2 = np.ones_like(aw_diff) * 1.5  # x-Werte für Gruppe 2
axs.plot(x2, aw_diff, color='lightcoral', label='AW', marker='s', linestyle = 'dotted', linewidth=0.25)

x3 = np.ones_like(ut_diff) * 1.78  # x-Werte für Gruppe 2
axs.plot(x3, ut_diff, color='lightgreen', label='UT', marker='D', linestyle = 'dotted', linewidth=0.25)


bars = ("Zape3", "Zape2", "Zape1","AW110A", "AW108", "AW107", "UT43.2", "UT30.3")
y_pos = np.arange(len(bars))
#bar1=axs[1].plot(y_pos, avg_read, zorder=3, linewidth=0.7)
axs.set_ylabel("% of reads filtered out")
#axs.set_title('b',loc='left', fontweight="bold")
axs.grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
axs.spines[['top','right']].set_visible(False)
#axs[1].bar_label(bar1, padding=3, color='black', fontsize=8)
axs.set_xticks([1.1, 1.22, 1.5, 1.78, 1.9])
axs.set_xticklabels(['','Zape', 'AW', 'UT',''])

plt.rcParams['figure.dpi']=500
plt.show()
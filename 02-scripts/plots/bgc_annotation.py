import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("/mnt/archgen/users/richtermeier/bachelor_thesis/05-results/ANNO_BGC_contig_taxclassification.tsv", sep='\t')
df = pd.DataFrame(data)

taxa = df['lineage'].str.extract(r'g_(.*?);')
taxa = taxa.iloc[:,0].value_counts()
taxa.rename('abundance', inplace=True)
taxa = pd.DataFrame(taxa)

fig, ax = plt.subplots()

ax.barh(taxa.index[0:20], taxa.iloc[0:20,0], color='cornflowerblue', zorder=3)
ax.set_xlabel("Amount of classified MAGs")
ax.grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax.spines[['top','right']].set_visible(False)
ax.set_title('b',loc='left', fontweight="bold")


plt.rcParams['figure.dpi']=300
plt.tight_layout()
plt.show()
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("/mnt/archgen/users/richtermeier/bachelor_thesis/05-results/ANNO_BGC_contig_taxclassification.tsv", sep='\t')
df = pd.DataFrame(data)

phylum = df['lineage'].str.extract(r'g_(.*?);')
phylum = phylum.iloc[:,0].value_counts()
phylum.rename('abundance', inplace=True)
phylum = pd.DataFrame(phylum)
fig, ax = plt.subplots()

ax.barh(phylum.index[0:20], phylum.iloc[0:20,0], color='cornflowerblue', zorder=3)
ax.set_xlabel("Amount of classified contigs with BGCs")
ax.grid(True, color='gray', linestyle=':', linewidth=0.5, zorder=0)
ax.spines[['top','right']].set_visible(False)


plt.rcParams['figure.dpi']=300
plt.tight_layout()
plt.show()
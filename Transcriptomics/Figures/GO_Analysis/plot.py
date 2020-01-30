import numpy as np
import pandas as pd
import hotspot.modules
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from tqdm import tqdm

# %% Load data

z_scores = pd.read_table(
    "../../CD4_w_protein/hotspot/hotspot_pairs2_z.txt.gz",
    index_col=0
)

hs_results = pd.read_table(
    "../../CD4_w_protein/hotspot/hotspot_hvg.txt",
    index_col=0
)

ens_map = {x: y for x, y in zip(hs_results.index, hs_results.Symbol)}

# %% Load Modules

modules = pd.read_table(
    "../../CD4_w_protein/hotspot/modules.txt",
    index_col=0
).Cluster

Z = pd.read_table(
    "../../CD4_w_protein/hotspot/linkage.txt",
    header=None
).values

# %% Modules to gene sets

gene_sets = {}
for i in modules.unique():
    if i == -1: continue
    genes = modules[modules == i].index
    genes = {ens_map[x].upper() for x in genes}
    gene_sets[i] = genes

all_genes = {ens_map[x].upper() for x in hs_results.index}

# %% Load GO sets

from gene_enrich import load_gene_set_gmt

# go_sets = load_gene_set_gmt("/data/yosef2/users/david.detomaso/Signatures/GO/GO_biological_process.gmt")
go_sets = load_gene_set_gmt("/data/yosef2/users/david.detomaso/Signatures/Enrichr/GO_Biological_Process_2015.txt")

# The Enrichr one is better!

# %% Overlap analysis
from  gene_enrich.analyses import gene_set_enrichment

results = {}
for k in tqdm(gene_sets):
    results[k] = gene_set_enrichment(
        gene_sets[k], all_genes, go_sets, correct_false_neg=True
    )

# %% Plot the top terms

from natsort import natsorted
top_terms = []
for m in natsorted(results)[::-1]:
    vals = results[m].iloc[0]
    top_terms.append(
        [str(m), vals.name, vals.pvalue, vals.genes]
    )
top_terms = pd.DataFrame(
    top_terms, columns=['Module', 'Name', 'FDR', 'Genes']
)

fig, ax1 = plt.subplots(figsize=(9, 5))
plt.barh(y=top_terms.Genes, width=np.log10(top_terms.FDR)*-1)
plt.yticks(size=8)
plt.xlabel('FDR ($-log_{10}$)')

ax2 = ax1.twinx()
labels = ['Module {}: {}'.format(x, y) for x, y in zip(top_terms.Module, top_terms.Name)]
labels = [x.split('(')[0] for x in labels]
plt.barh(y=labels, width=np.log10(top_terms.FDR)*-1)
plt.subplots_adjust(left=0.5, right=0.6, bottom=0.2)
plt.yticks(size=8)
plt.show()
# plt.savefig('GOTerms.svg')


# %% Why are 0 and 6 different?

g0 = modules.index[modules == 0]
g6 = modules.index[modules == 6]

z0 = z_scores.loc[g0, g0].values.ravel()
z6 = z_scores.loc[g6, g6].values.ravel()
z06 = z_scores.loc[g0, g6].values.ravel()


plt.figure()
kwargs = dict( 
    range=(-100, 100),
    alpha=0.5,
    bins=30,
)
plt.hist(z0, **kwargs, label='z0')
plt.hist(z6, **kwargs, label='z6')
plt.hist(z06, **kwargs, label='z06')
plt.legend()
plt.show()


zz06 = z_scores.loc[g0, g6]

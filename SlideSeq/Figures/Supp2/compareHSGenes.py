import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

plt.rcParams['svg.fonttype'] = 'none'

# %% print overlap stats

puck_analysis_dirs = {
    '9': '../../Puck_180819_9/',
    '10': '../../Puck_180819_10/',
    '11': '../../Puck_180819_11/',
    '12': '../../Puck_180819_12/',
}

def load_hs_genes(analysis_dir):
    results = pd.read_table(
        os.path.join(analysis_dir, "hotspot/hotspot.txt"),
        index_col=0
    )
    return results

hs_genes = {p: load_hs_genes(dd) for p, dd in puck_analysis_dirs.items()}

# %%

fout = open("overlap_stats.txt", "w")

fout.write(
    "{}\t{}\t{}\t{}\t{}\n".format(
        "Puck A", "Genes A", "Puck B", "Genes B", "Overlap"
    )
)

for a, b in itertools.combinations(hs_genes.keys(), 2):

    g1 = hs_genes[a].index[(hs_genes[a].FDR < .1)]
    g2 = hs_genes[b].index[(hs_genes[b].FDR < .1)]

    fout.write(
        "{}\t{}\t{}\t{}\t{}\n".format(
            a, g1.size, b, g2.size, (g1&g2).size
        )
    )

fout.close()

"""
Selects the genes using parameters from the hotspot output
"""

import os
import pandas as pd

hs_results_file = snakemake.input['hs_results']
genes_out = snakemake.output['genes']

os.makedirs(os.path.dirname(genes_out), exist_ok=True)

fdrThresh = snakemake.params['fdrThresh']

try:
    topN = int(snakemake.params['topN'])
except AttributeError:
    topN = None

try:
    cThresh = float(snakemake.params['cThresh'])
except AttributeError:
    cThresh = None

hs_results = pd.read_table(hs_results_file, index_col=0)


if topN is None:
    hs_genes = hs_results.index[hs_results.FDR < fdrThresh]
else:
    hs_genes = hs_results.sort_values('Z').tail(topN).index

if cThresh is not None:
    hs_genes = hs_results.index[hs_results.C > cThresh] & hs_genes

genes = list(hs_genes)

with open(genes_out, 'w') as fout:
    for g in genes:
        fout.write(g + "\n")

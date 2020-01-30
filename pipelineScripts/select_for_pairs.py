import os
import pandas as pd

gene_effect_file = snakemake.input['ge_file']
gene_effects = pd.read_table(
    gene_effect_file, index_col=0
)
genes_out = snakemake.output['genes']

os.makedirs(os.path.dirname(genes_out), exist_ok=True)


ge_size = gene_effects.iloc[:, 1:].abs().max(axis=1)
ge_size = ge_size.loc[ge_size > 0]

genes = ge_size.index.tolist()

with open(genes_out, 'w') as fout:
    for g in genes:
        fout.write(g + "\n")

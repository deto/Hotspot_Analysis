import os
import loompy

loom_file = snakemake.input['loom']
genes_out = snakemake.output['genes']

os.makedirs(os.path.dirname(genes_out), exist_ok=True)

N_CELLS = snakemake.params['N_CELLS']

with loompy.connect(loom_file, 'r') as ds:
    ens_id = ds.ra['EnsID'][:]
    counts = ds[:, :]

cell_counts = (counts > 0).sum(axis=1)

gene_passes = cell_counts >= N_CELLS
genes = list(ens_id[gene_passes])

with open(genes_out, 'w') as fout:
    for g in genes:
        fout.write(g + "\n")

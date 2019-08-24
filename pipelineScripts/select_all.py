import os
import loompy

loom_file = snakemake.input['loom']
genes_out = snakemake.output['genes']

os.makedirs(os.path.dirname(genes_out), exist_ok=True)

with loompy.connect(loom_file, 'r') as ds:
    ens_id = ds.ra['EnsID'][:]

genes = list(ens_id)

with open(genes_out, 'w') as fout:
    for g in genes:
        fout.write(g + "\n")

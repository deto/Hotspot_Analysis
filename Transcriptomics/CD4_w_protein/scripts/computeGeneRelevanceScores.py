import os
import gene_enrich
import pandas as pd


sig_file = snakemake.input['signatures']
out_file = snakemake.output['out']

out_dir = os.path.dirname(out_file)
os.makedirs(out_dir, exist_ok=True)

signatures = gene_enrich.load_gene_set_gmt(sig_file)

# Count all genes

from collections import Counter

all_genes = []
for x in signatures:
    all_genes.extend(x.genes)

gene_counts = {gene: count for gene, count in Counter(all_genes).items()}
gene_counts = (pd.Series(gene_counts)
               .sort_values(ascending=False)
               .rename('Counts')
               )

gene_counts.to_frame().to_csv(out_file, sep="\t")

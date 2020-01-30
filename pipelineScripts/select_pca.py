import pandas as pd

comp_file = snakemake.input['components']
genes_out = snakemake.output['genes']
gene_info_out = snakemake.output['geneInfo']

comps = pd.read_table(comp_file, index_col=0)

gene_scores = comps.abs().sum(axis=1)

genes = gene_scores.sort_values(ascending=False).index[0:500].tolist()

with open(genes_out, 'w') as fout:
    for g in genes:
        fout.write(g + "\n")

gene_scores.name = "Score"
gene_scores.to_frame().to_csv(gene_info_out, sep="\t")

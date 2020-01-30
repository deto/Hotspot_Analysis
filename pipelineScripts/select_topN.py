import pandas as pd

gene_info_file = snakemake.input['geneInfo']

N = snakemake.params['N']
variable = snakemake.params['var']
ascending = snakemake.params['ascending']

genes_out = snakemake.output['genes']

gene_info = pd.read_table(gene_info_file, index_col=0)

gene_info = gene_info.sort_values(variable, ascending=ascending)

genes = gene_info.index[0:N]

with open(genes_out, 'w') as fout:
    for g in genes:
        fout.write(g + "\n")


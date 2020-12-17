import loompy
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import time

import sys


# loom_file = snakemake.input['loom']
# genes_file = snakemake.input['genes']
# out_file_scores = snakemake.output['pair_scores']

if __name__ == '__main__':

    loom_file = sys.argv[1]
    genes_file = sys.argv[2]
    out_file_scores = sys.argv[3]

    print(loom_file)
    print(genes_file)
    print(out_file_scores)

    loom_file = os.path.abspath(loom_file)
    genes_file = os.path.abspath(genes_file)

    with loompy.connect(loom_file, 'r') as ds:
        barcodes = ds.ca['Barcode'][:]
        counts = ds[:, :]
        gene_info = ds.ra['EnsID', 'Symbol']
        num_umi = ds.ca['NumUmi'][:]


    valid_genes = pd.read_table(genes_file, header=None)[0].values

    # Have to do this because data_slideseq makes it a numpy array
    gene_info = pd.DataFrame(
        gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
    counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)

    from arboreto.algo import grnboost2, genie3

    log_scaled_counts = (
        np.log(counts.divide(counts.sum(axis=0), axis=1) * 10000 + 1)
    )

    log_scaled_counts = log_scaled_counts.loc[valid_genes, :]

    # About 14 minutes for 5000 genes x 3000 cells
    old_dir = os.getcwd()
    os.makedirs("/scratch/david.detomaso/temp", exist_ok=True) # Need this or else the workers time out
    os.chdir("/scratch/david.detomaso/temp")

    a = time.time()
    network = grnboost2(log_scaled_counts.T)
    b = time.time()

    print(b-a)

    os.chdir(old_dir)

    # Need to convert to long

    net_wide = network.pivot(index='TF', columns='target', values='importance')
    z = net_wide.fillna(0)
    z = z + z.T

    z.to_csv(out_file_scores, sep="\t", compression="gzip")

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

import NaiveDE, SpatialDE


# Load some data

latent_file = "/data/yosef2/users/david.detomaso/VISION/SlideSeq/Puck_180819_12/BeadLocationsForR_filtered.csv"
latent = pd.read_csv(latent_file, index_col=0)


genes_file = "/data/yosef2/users/david.detomaso/VISION/SlideSeq/Puck_180819_12/genes.txt"
genes = open(genes_file).readlines()
genes = [x.strip() for x in genes]

data_file = "/data/yosef2/users/david.detomaso/VISION/SlideSeq/Puck_180819_12/data.npz"
data = np.load(data_file)

counts = csc_matrix(
    (data['data'], data['indices'], data['indptr']),
    shape=data['shape']
)

counts = counts.T.toarray()

num_umi = counts.sum(axis=0)

# Do some filtering
valid_genes = (counts > 0).sum(axis=1) > 10
counts = counts[valid_genes, :]
genes = [x for x, y in zip(genes, valid_genes) if y]

counts = counts.T

counts = pd.DataFrame(counts, columns=genes)

# prepare 'sample_info' like they want it
sample_info = latent.copy()
sample_info['total_counts'] = num_umi

# Some normalization

import time

a = time.time()
norm_expr = NaiveDE.stabilize(counts.T).T
b = time.time()
print(b-a)  # 13 seconds for 30k cells, 11k genes

a = time.time()
resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T
b = time.time()
print(b-a) # 4 seconds for 30k cells, 11k, genes


sample_resid_expr = resid_expr.sample(n=20, axis=1, random_state=1)
X = sample_info[['xcoord', 'ycoord']]

a = time.time()
results = SpatialDE.run(X, sample_resid_expr)
b = time.time()
print(b-a) # 10460 seconds work 30k cells, 10 genes
print(b-a) # 10356 seconds work 30k cells, 20 genes

# Now do it for everything
a = time.time()
results = SpatialDE.run(X, resid_expr)
b = time.time()
print(b-a) # 26347 seconds for 30k cells, 11,051 genes
results = results.set_index('g')
results.to_csv("spatial_de_results.csv")


# Now try the histology
# qvalues are...weird for some reason
# Seem to be MORE inflated than pvalues
# Just take an extreme value to get around 500 genes
sign_results = results.reset_index().query('qval < 0.0001')

sign_results['l'].value_counts()

#  388.377380     255
#  118.866849     216
#  1268.957582     72
#  36.380409       19
#  3.407856         7
#  4146.104863      2
#  Name: l, dtype: int64


# Results with really high l (1268) seem to be garbage.  Choose 350 instead 

a = time.time()
histology_results, patterns = SpatialDE.aeh.spatial_patterns(
    X, resid_expr, sign_results, C=3, l=350, verbosity=1)
b = time.time()
print(b-a)  # I killed this off after 1.5 days.
# It was only on ITER 1 out of maybe 8 or so

histology_results.to_csv("spatial_de_histology.csv")
patterns.to_csv("spatial_de_patterns.csv")

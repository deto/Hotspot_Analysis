import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib
import loompy
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MAX_COMPONENTS = 100


loom_file = snakemake.input['loom']
pca_out = snakemake.output['latent']

out_dir = os.path.dirname(pca_out)


with loompy.connect(loom_file, 'r') as ds:
    expression = ds.layers['scaled'][:]
    ens_ids = ds.ra['EnsID'][:]
    barcodes = ds.ca['Barcode'][:]

expression = pd.DataFrame(expression, index=ens_ids, columns=barcodes)
expression = np.log2(expression + 1)

# Gene filtering
genes = snakemake.input['genes']
genes = pd.read_table(genes, header=None).iloc[:, 0].tolist()
expression = expression.loc[genes]


model = PCA(n_components=MAX_COMPONENTS, random_state=0)
model.fit(expression.values.T)


# Now, look at null data
null = expression.values.copy()
for i in range(null.shape[0]):
    np.random.shuffle(null[i, :])

null_model = PCA(n_components=MAX_COMPONENTS, random_state=0)
null_model.fit(null.T)


# First value less than = num PCs
# E.g., if index 21 is first value less than null, then 0-20
# are greater than null, and num_pcs = 21
# num_pcs = np.nonzero(
#     (model.explained_variance_ < null_model.explained_variance_)
# )[0][0]

num_pcs = 20

# Diagnostic plot
xmin = 0
xmax = min(num_pcs*2, 200)
ymin = 0
ymax = max(model.explained_variance_.max(),
           null_model.explained_variance_.max())*1.2

plt.figure()
plt.plot(np.arange(len(model.explained_variance_))+1,
         model.explained_variance_, 'o-', label='Expression')
plt.plot(np.arange(len(null_model.explained_variance_))+1,
         null_model.explained_variance_, 'o-', label='Null')

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

plt.vlines(num_pcs, ymin, ymax, linestyle='dashed', color='red')
plt.legend(loc='best')
plt.title('Significant PCs: {}'.format(num_pcs))
plt.savefig(os.path.join(out_dir, 'PermutationPA.png'))

pcs = model.fit_transform(expression.values.T)
pcs = pcs[:, 0:num_pcs]

pcs = pd.DataFrame(pcs, index=expression.columns)

pcs.to_csv(pca_out, sep='\t')

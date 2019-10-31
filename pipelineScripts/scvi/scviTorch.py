import os

import numpy as np
import loompy
from sklearn.manifold import TSNE

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import torch
from scvi.dataset import CortexDataset, RetinaDataset
from scvi.dataset.dataset import GeneExpressionDataset

from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer

import pandas as pd

#gene_dataset = CortexDataset(save_path="cortex")

loom_file = snakemake.input['loom']
genes = snakemake.input['genes']

latent_file = snakemake.output['latent']
model_file = snakemake.output['model']

genes = pd.read_table(genes, header=None).iloc[:, 0].tolist()

out_dir = os.path.dirname(latent_file)

os.makedirs(out_dir, exist_ok=True)

with loompy.connect(loom_file, 'r') as ds:
    counts = ds.layers[''][:]
    ens_ids = ds.ra['EnsID'][:]
    barcodes = ds.ca['Barcode'][:]

counts = pd.DataFrame(counts, index=ens_ids, columns=barcodes)

# Do some filtering
counts = counts.loc[genes]

counts = counts.astype("int32").T

batch = np.zeros(counts.shape[0])

# Save batch encoding for later loading

try:
    COMPONENTS = int(snakemake.params['components'])
except AttributeError:
    COMPONENTS = 10  # Latent component count

try:
    LAYERS = int(snakemake.params['layers'])
except AttributeError:
    LAYERS = 1  # number of hidden layers

RECONSTRUCTION_LOSS = "nb"


batch = batch.reshape((-1, 1)).astype('int64')
cvals = counts.values.astype('int64')
zz = GeneExpressionDataset.get_attributes_from_matrix(
    cvals, batch_indices=batch)

# zz[0]: int64, ndarray, genes x cells
# zz[1]: float32, ndarray, cells x 1
# zz[2]: float32, ndarray, cells x 1
# zz[3]: int64, ndarray, cells x 1

dataset = GeneExpressionDataset(*zz, gene_names=counts.columns)

n_epochs = 400
lr = 1e-3
use_batches = True
use_cuda = torch.cuda.is_available()

torch.set_num_threads(20)
# However, need to set MKL_NUM_THREADS too

print("Running scVI with {} components, on {} batches and {} genes".format(
    COMPONENTS, dataset.n_batches, cvals.shape[1])
)

vae = VAE(dataset.nb_genes, n_batch=dataset.n_batches * use_batches,
          n_latent=COMPONENTS, n_layers=LAYERS,
          reconstruction_loss=RECONSTRUCTION_LOSS)
trainer = UnsupervisedTrainer(vae,
                              dataset,
                              train_size=0.75,
                              use_cuda=use_cuda,
                              frequency=5)
trainer.train(n_epochs=n_epochs, lr=lr)


# Plot training result

plt.figure()
ll_train_set = trainer.history["ll_train_set"]
ll_test_set = trainer.history["ll_test_set"]
x = np.linspace(0, n_epochs, (len(ll_train_set)))
plt.plot(x, ll_train_set)
plt.plot(x, ll_test_set)
ymax = np.percentile(ll_test_set, 95)
ymin = np.min(ll_train_set) - .5*(ymax - np.min(ll_train_set))
plt.ylim(ymin, ymax)
plt.savefig(os.path.join(out_dir, 'training.png'))

tData = pd.DataFrame({
    'Training_Loss': ll_train_set,
    'Test_Loss': ll_test_set,
})

tData.to_csv(os.path.join(out_dir, 'training.txt'), sep="\t")


# Get latent space
latent = trainer.get_all_latent_and_imputed_values()["latent"]

# Save Results
latent = pd.DataFrame(latent, index=counts.index)
latent.to_csv(latent_file, sep="\t", compression="gzip")

torch.save(trainer.model, model_file)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loompy
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import NearestNeighbors
from scipy.stats import multivariate_normal, poisson, norm
from sklearn.decomposition import PCA
from tqdm import tqdm
from sklearn.manifold import TSNE


def compute_correlation_matrix(coordinates, bandwidth_k=30):

    nn = NearestNeighbors(n_neighbors=bandwidth_k)
    nn.fit(coordinates)
    dist, ind = nn.kneighbors()

    bandwidth = np.median(dist[:, -1])

    dist_mat = squareform(pdist(coordinates))

    corr_mat = np.exp(-1 * (dist_mat/bandwidth)**2)

    return corr_mat


def compute_correlation_matrix_rank(coordinates, bandwidth_k=30):
    """
    Uses the ranking of cells distances to each other to determine
    their correlation

    This doesn't work because the resulting correlation matrix
    is not positive semidefinite
    """

    dist_mat = squareform(pdist(coordinates))
    dist_ranks = np.argsort(np.argsort(dist_mat, axis=1), axis=1)

    dist_ranks = np.minimum(dist_ranks, dist_ranks.T)

    corr_mat = np.exp(-1 * dist_ranks / bandwidth_k)

    return corr_mat


def simulate_autocorrelated_counts(corr_mat, gene_means):
    """
    gene_means are in counts/10k scale
    """
    N_CELLS = corr_mat.shape[0]

    gene_mean = np.random.choice(gene_means)
    bcv = np.random.uniform(.5, 2, 1)[0]
    sd = bcv * gene_mean

    # convert to log
    gene_mean_log = np.log10(gene_mean)
    gene_sd_log = np.log10(gene_mean+sd) - gene_mean_log

    per_cell_mu = np.ones(N_CELLS) * gene_mean_log

    cell_means = multivariate_normal(
        mean=per_cell_mu, cov=corr_mat * (gene_sd_log**2), allow_singular=True
    ).rvs(1)

    cell_counts = poisson.rvs(
        10**cell_means, size=N_CELLS
    )

    return cell_counts


def simulate_autocorrelated_counts_from_component(component, gene_means):
    """
    Given a perturbation component, simulate the counts for one random gene
    """
    N_CELLS = corr_mat.shape[0]

    gene_mean = np.random.choice(gene_means)
    bcv = np.random.uniform(.5, 2, 1)[0]
    sd = bcv * gene_mean

    # convert to log
    gene_mean_log = np.log10(gene_mean)
    gene_sd_log = np.log10(gene_mean+sd) - gene_mean_log

    per_cell_mu = np.ones(N_CELLS) * gene_mean_log

    cell_means = component * (gene_sd_log**2) + per_cell_mu

    cell_counts = poisson.rvs(
        10**cell_means, size=N_CELLS
    )

    return cell_counts


def simulate_autocorrelated_components(corr_mat, K):
    """
    Given the corr_mat, simulate components of variation multivariate

    Return as a cells x K matrix
    """

    n_cells = corr_mat.shape[0]
    per_cell_mu = np.zeros(n_cells)

    components = multivariate_normal(
        mean=per_cell_mu, cov=corr_mat, allow_singular=True
    ).rvs(K).T

    return components


def simulate_null_counts(corr_mat, gene_means):
    """
    gene_means are in counts/10k scale
    """
    N_CELLS = corr_mat.shape[0]

    gene_mean = np.random.choice(gene_means)
    bcv = np.random.uniform(.5, 2, 1)[0]
    sd = bcv * gene_mean

    # convert to log
    gene_mean_log = np.log10(gene_mean)
    gene_sd_log = np.log10(gene_mean+sd) - gene_mean_log

    cell_means = gene_mean_log

    cell_counts = poisson.rvs(
        10**cell_means, size=N_CELLS
    )

    return cell_counts


def load_reference_info(reference_dataset_expression_loom):

    with loompy.connect(reference_dataset_expression_loom, 'r') as ds:
        counts = ds[:, :]

    counts = counts[counts.sum(axis=1) > 10]
    counts_scaled = counts / counts.sum(axis=0, keepdims=True) * 10000
    gene_means = counts_scaled.mean(axis=1)

    model = PCA(n_components=20)
    pca = model.fit_transform(np.log2(counts_scaled.T + 1))

    return gene_means, pca

# %%


# reference_dataset_expression_loom = "../10x_PBMC_w_proteins/cd4/data.loom"
reference_dataset_expression_loom = snakemake.input['ref_loom']
out_file = snakemake.output['out_loom']
out_pca = snakemake.output['out_pca']
out_tsne = snakemake.output['out_tsne']
out_components = snakemake.output['out_components']

N_components = 5

gene_means, pca = load_reference_info(reference_dataset_expression_loom)

model = TSNE(n_components=2)
tsne = model.fit_transform(pca)

corr_mat = compute_correlation_matrix(tsne, 100)

components = simulate_autocorrelated_components(corr_mat, K=N_components)

N_pos_genes = 1000
N_null_genes = 9000
N_cells = corr_mat.shape[0]

out_pos_counts = np.zeros((N_pos_genes, N_cells))
out_null_counts = np.zeros((N_null_genes, N_cells))

gene_component_assignments = np.zeros(N_pos_genes + N_null_genes)

for i in tqdm(range(N_pos_genes)):
    component_i = np.random.choice(N_components)
    out_pos_counts[i] = simulate_autocorrelated_counts_from_component(
        components[:, component_i], gene_means)
    gene_component_assignments[i] = component_i

for i in tqdm(range(N_null_genes)):
    out_null_counts[i] = simulate_null_counts(corr_mat, gene_means)
    gene_component_assignments[i] = -1

out_symbols = ['PosGene{}'.format(i+1) for i in range(N_pos_genes)] + \
              ['NullGene{}'.format(i+1) for i in range(N_null_genes)]

out_ensid = out_symbols

out_barcodes = ['Cell{}'.format(i+1) for i in range(N_cells)]

out_counts = np.concatenate((out_pos_counts, out_null_counts), axis=0)
out_num_umi = out_counts.sum(axis=0)

out_scaled = out_counts / out_counts.sum(axis=0, keepdims=True) * 10000

out_pos_gene = np.array(
    [True]*N_pos_genes + [False]*N_null_genes
)

row_attrs = {
    "Symbol": np.array(out_symbols),
    "EnsID": np.array(out_ensid),
    "PosGene": out_pos_gene,
    "GeneComponent": gene_component_assignments,
}

col_attrs = {
    "Barcode": np.array(out_barcodes),
    "NumUmi": out_num_umi,
}

layers = {
    '': out_counts,
    'scaled': out_scaled
}

loompy.create(out_file, layers, row_attrs, col_attrs)


pca = pd.DataFrame(
    pca,
    index=out_barcodes,
    columns=["PC{}".format(i + 1) for i in range(pca.shape[1])],
)

pca.to_csv(out_pca, sep="\t")

tsne = pd.DataFrame(
    tsne,
    index=out_barcodes,
    columns=["TSNE{}".format(i + 1) for i in range(tsne.shape[1])],
)

tsne.to_csv(out_tsne, sep="\t")

components = pd.DataFrame(
    components,
    index=out_barcodes,
    columns=["Component{}".format(i) for i in range(components.shape[1])],
)

components.to_csv(out_components, sep="\t")

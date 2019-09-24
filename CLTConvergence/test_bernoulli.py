"""
Testing the limits on the convergence of the bernoulli model to the
assumptions put forward by the Central Limit Theorem

This is an easier case than DANB because there is only one
distributional parameter to scan.
"""

from hotspot.sim_data import (
    sim_latent, sim_counts_bernoulli, sim_umi_counts,
    generate_permutation_null
)
from hotspot import Hotspot
from tqdm import tqdm
import numpy as np
import pandas as pd


N_CELLS = 30000
N_NEIGHBORS = 30
N_LATENT_DIM = 20
N_UMI = 50

latent = sim_latent(N_CELLS, N_LATENT_DIM)
num_umi = np.ones(N_CELLS) * N_UMI

num_umi = pd.Series(num_umi)
latent = pd.DataFrame(latent)

N_REPS = 10000

gene_p_trials = [.003, .01, .03, .1, .3, 1, 3, 10, 30]

trial_results = []
trial_Gs = []

for gene_p in gene_p_trials:

    detect_p = 1 - (1-gene_p/10000)**N_UMI
    N_DETECTS = round(detect_p * N_CELLS)

    if N_DETECTS < 1:
        continue

    cc = np.array([1]*N_DETECTS + [0]*(N_CELLS-N_DETECTS))
    counts = generate_permutation_null(cc, N_REPS)

    counts = pd.DataFrame(counts)

    hs = Hotspot(counts, latent, num_umi)
    hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
    results = hs.compute_hotspot(model='bernoulli', centered=True, jobs=20)

    fp = (results.Pval < .05).sum()
    fdr_rate = (results.FDR < .05).sum()

    trial_results.append([gene_p, fp, fdr_rate, results.stdG[0]])
    trial_Gs.append(results.G)

trial_results = pd.DataFrame(
    trial_results, columns=['gene_p', 'fp', 'fdr_rate', 'stdG']
)

plt.figure()
plt.hist(trial_Gs[2], 50)
plt.show()

gene_p = .7
detect_p = 1 - (1-gene_p/10000)**N_UMI
N_DETECTS = round(detect_p * N_CELLS)
detect_p = N_DETECTS/N_CELLS

# Approx number of coincident detections:
from scipy.stats import binom
cd = N_CELLS * detect_p * binom.mean(n=N_NEIGHBORS, p=(N_DETECTS-1)/N_CELLS)
cd = N_CELLS * detect_p * binom.mean(n=N_NEIGHBORS, p=detect_p - 1/N_CELLS)
cd = N_CELLS * detect_p * N_NEIGHBORS * (detect_p - 1/N_CELLS)
cd = N_CELLS * N_NEIGHBORS * detect_p**2 - N_NEIGHBORS * detect_p

# approximately
cd = N_CELLS * N_NEIGHBORS * detect_p**2

# suggests for cd = 10:
detect_p_thresh = (10 / (N_CELLS * N_NEIGHBORS))**0.5
gene_p_thresh = (1-(1-detect_p_thresh)**(1/N_UMI))*10000

n_detect_thresh = detect_p_thresh * N_CELLS
n_detect_thresh = (10 * N_CELLS / N_NEIGHBORS)**.5

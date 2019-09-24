"""
For some reason, the DANB model doesn't result in correct variance estimates
for the autocorrelation statistic.

Here...exploring why that happens.

- Actually, as far as the centered model goes, it doesn't appear to happen with simulated data.
- Permuted data is ok too (with centered model)

What about non-centered?
    - we appear to be off by exactly 0.5?  That's....odd
"""

from hotspot.sim_data import (
    sim_latent, sim_umi_counts, sim_counts_danb,
    generate_permutation_null
)
from hotspot import Hotspot
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hotspot import danb_model
from importlib import reload
from hotspot import local_stats
from hotspot.sim_data import sim_counts_bernoulli
from numba import jit


N_CELLS = 5000
N_NEIGHBORS = 30
N_LATENT_DIM = 20
N_UMI = 2000

latent = sim_latent(N_CELLS, N_LATENT_DIM)
num_umi = np.ones(N_CELLS) * N_UMI

num_umi = pd.Series(num_umi)
latent = pd.DataFrame(latent)

with_neg_size = []
for i in tqdm(range(20)):
    counts = sim_counts_danb(N_CELLS, 1000, num_umi.values, mean=10, size=10000)
    counts = pd.DataFrame(counts)

    hs = Hotspot(counts, latent, num_umi)
    hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
    results = hs.compute_hotspot(model='danb', centered=True, jobs=20)

    with_neg_size.append(results.G.std())

without_neg_size = []
for i in tqdm(range(20)):
    counts = sim_counts_danb(N_CELLS, 1000, num_umi.values, mean=10, size=10000)
    counts = pd.DataFrame(counts)

    hs = Hotspot(counts, latent, num_umi)
    hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
    results = hs.compute_hotspot(model='danb', centered=True, jobs=20)

    without_neg_size.append(results.G.std())

counts = sim_counts_danb(N_CELLS, 1000, num_umi.values, mean=10, size=.2)
counts = pd.DataFrame(counts)

hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=True, jobs=20)
print(results.G.std(), results.stdG[0])


plt.figure()
plt.hist(with_neg_size, bins=20, range=(430, 470), alpha=0.5)
plt.hist(without_neg_size, bins=20, range=(430, 470), alpha=0.5)
plt.show()


plt.figure()
plt.hist(results.G, bins=100)
plt.show()


plt.figure()
plt.plot(
    results.G,
    counts.var(axis=1) / counts.mean(axis=1),
    'o', ms=2
)
plt.show()



sizes = []
mus = []
for i in tqdm(range(counts.shape[0])):
    mu, var, x2 = danb_model.fit_gene_model(counts.iloc[i, :].values, num_umi.values)
    size = mu[0] / (var[0]/mu[0] - 1)
    sizes.append(size)
    mus.append(mu[0])

plt.figure()
#plt.plot(results.G, sizes, 'o', ms=2)
plt.plot(results.G, mus, 'o', ms=2)
plt.show()

plt.figure()
plt.plot(mus, sizes, 'o', ms=2)
plt.show()

reload(danb_model)

# Now, what about permuted data?
counts1 = sim_counts_danb(N_CELLS, 1, num_umi.values, mean=1, size=.5).ravel()
counts = generate_permutation_null(counts1, N_REPS=1000)
counts = pd.DataFrame(counts)

sizes = []
mus = []
for i in tqdm(range(counts.shape[0])):
    mu, var, x2 = danb_model.fit_gene_model(counts.iloc[i, :].values, num_umi.values)
    size = mu[0] / (var[0]/mu[0] - 1)
    sizes.append(size)
    mus.append(mu[0])

hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=True, jobs=20)

# What about non-centered model?
counts1 = sim_counts_danb(N_CELLS, 1, num_umi.values, mean=10, size=1000).ravel()
counts = generate_permutation_null(counts1, N_REPS=1000)
counts = pd.DataFrame(counts)

hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=False, jobs=20)
print(results.Z.std())

counts = sim_counts_danb(N_CELLS, 1000, num_umi.values, mean=10, size=1000)
counts = pd.DataFrame(counts)
hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=False, jobs=20)
print(results.Z.std())

# Bernoulli has the same issue when it's non-centered
counts1 = sim_counts_bernoulli(N_CELLS, num_umi.values, gene_p=.1)
counts = generate_permutation_null(counts1, N_REPS=1000)
counts = pd.DataFrame(counts)
hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='bernoulli', centered=False, jobs=20)
print(results.Z.std())

# is the computation correct?
from hotspot import local_stats
reload(local_stats)
counts = sim_counts_danb(N_CELLS, 1000, num_umi.values, mean=10, size=1000)
counts = pd.DataFrame(counts)

hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=False, jobs=20)

mu, var, x2 = danb_model.fit_gene_model(counts.iloc[33, :].values, num_umi.values)

eg, eg2 = local_stats.compute_moments_weights_slow(
    mu, x2, hs.neighbors.values, hs.weights.values)

egF, eg2F = local_stats.compute_moments_weights(
    mu, x2, hs.neighbors.values, hs.weights.values)

stdG = (eg2 - eg**2)**.5
print(stdG)
stdGF = (eg2F - egF**2)**.5
print(stdGF)

# One observation - they are slightly different
# However, difference is very small 1 part in 10,000 roughly
# Too small to explain what we are seeing

# Instead, the real issue seems to be that the estimate of the MEAN is off

plt.figure()
plt.hist(results.Pval, 30)
plt.show()

C = (results.G - results.EG)


eg_original = N_CELLS*N_NEIGHBORS*mu[0]**2
print(eg_original)
eg_dep = N_CELLS*N_NEIGHBORS*(mu[0]**2 - var[0] / (N_CELLS - 1))
print(eg_dep)


xx = (results.G**2) - 2*results.G*results.EG + results.EG**2
xx.mean()**.5

a = (results.EG**2).mean()
b = (results.G*results.EG).mean()

results['EG2'] = results.stdG**2 + results.EG**2
xx = results.EG2 - 2*results.EG*results.EG + results.EG**2
xx.mean()**.5


np.log10((results.G*results.EG - results.EG**2).mean())



results.G.var()

(results.stdG**2).mean()


(results.EG2 - results.EG**2).mean()


(results.G - results.EG).var()


# Correlation between G and EG break this!!
cor = np.corrcoef(results.G, results.EG)[0, 1]
cov = cor * results.G.std() * results.EG.std()


(results.EG2*5000/4999 - results.EG**2)**.5

# But what about if we are just doing permutations?
# Under a permutation, E[G] is fixed...


counts1 = sim_counts_danb(N_CELLS, 1, num_umi.values, mean=10, size=1000).ravel()
counts = generate_permutation_null(counts1, N_REPS=1000)
counts = pd.DataFrame(counts)

mu, var, x2 = danb_model.fit_gene_model(counts.iloc[0, :].values, num_umi.values)

hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=False, jobs=20)
print(results.Z.std())

# I hate this so much



@jit(nopython=True)
def order3(x, ITER):
    tot = 0
    
    for i in range(ITER):
        nodes = np.random.choice(x.size, size=3, replace=False)
        tot += x[nodes[0]]**2*x[nodes[1]]*x[nodes[2]]

    return tot / ITER


x = counts.values[0]

o3s = []
for i in tqdm(range(100)):
    o3 = order3(x, 500000)
    o3s.append(o3)

i = np.arange(N_CELLS)

@jit(nopython=True)
def test(ix, x):

    N = ix.shape[1]

    tot = 0
    for i in range(N):
        tot += x[ix[0, i]] ** 2 * x[ix[1, i]] * x[ix[2, i]]

    return tot / N

ix = generate_permutation_null(i, 3)
test(ix, x)


# The null is wrong!  It's not over the space of graphs.
# Rather, it's over the space of K-nearest-neighbor graphs!

# Try this as the null:

@jit(nopython=True)
def compute_moments_weights_slow2(mu, x2, neighbors, weights, c):
    """
    This version exaustively iterates over all |E|^2 terms
    to compute the expected moments exactly.  Used to test
    the more optimized formulations that follow

    Here we assume that the FIRST neighbor in a pair is NOT random
    """

    N = neighbors.shape[0]
    K = neighbors.shape[1]

    # Calculate E[G]
    EG = 0
    for i in range(N):
        for k in range(K):
            j = neighbors[i, k]
            wij = weights[i, k]

            EG += wij*c[i]*mu[j]

    # Calculate E[G^2]
    count1 = 0
    count2 = 0
    count3 = 0
    EG2 = 0
    for i in range(N):

        EG2_i = 0

        for k in range(K):
            j = neighbors[i, k]
            wij = weights[i, k]

            for x in range(N):
                for z in range(K):

                    # Two edges (i, j) and (x, y)
                    
                    y = neighbors[x, z]
                    wxy = weights[x, z]

                    s = wij*wxy
                    if s == 0:
                        continue

                    if i == x:

                        if j == y:

                            t1 = c[i]**2 * x2[j]
                            count1 += 1

                        else:

                            t1 = c[i]**2 * mu[j] * mu[y]
                            count2 += 1

                    else:

                        t1 = c[i]*mu[j]*c[x]*mu[y]
                        count3 += 1

                    EG2_i += s * t1

        EG2 += EG2_i

    return EG, EG2, count1, count2, count3


neighbors = hs.neighbors.values
weights_x = np.ones_like(hs.weights.values)
c = counts.iloc[0, :].values

mu, var, x2 = danb_model.fit_gene_model(c, num_umi.values)

EG, EG2, c1, c2, c3 = compute_moments_weights_slow2(mu, x2, neighbors, weights_x, c)

# Does this describe a permutation over graphs?
# Try it:


@jit(nopython=True)
def gen_neighbors(N_CELLS, N_NEIGHBORS):

    nn = np.zeros((N_CELLS, N_NEIGHBORS), dtype=np.int64)

    for i in range(N_CELLS):
        vals = np.random.choice(N_CELLS, size=N_NEIGHBORS, replace=False)
        nn[i, :] = vals

    return nn

weights = np.ones((N_CELLS, N_NEIGHBORS))


vals = []
for i in tqdm(range(10000)):
    nn = gen_neighbors(N_CELLS, N_NEIGHBORS)
    G = local_stats.local_cov_weights(c, nn, weights)
    vals.append(G)

vals = np.array(vals)
vals.std()

(EG2 - EG**2)**.5


@jit(nopython=True)
def calc_degree_dist(neighbors, weights):

    N = neighbors.shape[0]
    K = neighbors.shape[1]

    D = np.zeros(N)

    for i in range(N):
        for k in range(K):

            j = neighbors[i, k]
            wij = weights[i, k]

            D[i] += wij
            D[j] += wij

    return D

D1 = calc_degree_dist(hs.neighbors.values, hs.weights.values)

nn = gen_neighbors(N_CELLS, N_NEIGHBORS)
weights = np.ones((N_CELLS, N_NEIGHBORS))
D2 = calc_degree_dist(nn, weights)

fig, axs = plt.subplots(1, 2)
plt.sca(axs[0])
plt.hist(D1, 30)
plt.title('Simulated Latent Space')
plt.sca(axs[1])
plt.hist(D2, 30)
plt.title('Randomly chosen neighbors')
plt.show()

EG, EG2, c1, c2, c3 = compute_moments_weights_slow2(
    mu, x2, hs.neighbors.values, hs.weights.values, c)


# Does this estimate the variance with the given degree distribution?
counts1 = sim_counts_danb(N_CELLS, 1, num_umi.values, mean=10, size=1000).ravel()
counts = generate_permutation_null(counts1, N_REPS=1000)
counts = pd.DataFrame(counts)

mu, var, x2 = danb_model.fit_gene_model(counts.iloc[0, :].values, num_umi.values)

hs_n = Hotspot(counts, latent, num_umi)
hs_n.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
hs_n.neighbors = pd.DataFrame(gen_neighbors(N_CELLS, N_NEIGHBORS))
hs_n.weights = pd.DataFrame(np.ones((N_CELLS, N_NEIGHBORS)))
results = hs_n.compute_hotspot(model='danb', centered=False, jobs=20)
print(results.Z.std())


def generate_resampled_null(observed_counts, N_REPS):

    out = []
    for i in range(N_REPS):
        out.append(
            np.random.choice(
                observed_counts, size=len(observed_counts), replace=True
            )
        )

    return np.vstack(out)


counts1 = sim_counts_danb(N_CELLS, 1, num_umi.values, mean=10, size=1000).ravel()
counts = generate_resampled_null(counts1, N_REPS=1000)
counts = pd.DataFrame(counts)

mu, var, x2 = danb_model.fit_gene_model(counts.iloc[0, :].values, num_umi.values)

hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=False, jobs=20)
print(results.Z.std())

# Try again with the permutations but compute EG based on degree distribution

counts1 = sim_counts_danb(N_CELLS, 1, num_umi.values, mean=10, size=1000).ravel()
counts = generate_permutation_null(counts1, N_REPS=1000)
counts = pd.DataFrame(counts)

mu, var, x2 = danb_model.fit_gene_model(counts.iloc[0, :].values, num_umi.values)

hs = Hotspot(counts, latent, num_umi)
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)
results = hs.compute_hotspot(model='danb', centered=False, jobs=20)
print(results.Z.std())

@jit(nopython=True)
def new_eg(x, mu, neighbors, weights):

    tot = 0
    N = neighbors.shape[0]
    K = neighbors.shape[1]

    for i in range(N):
        for k in range(K):

            j = neighbors[i, k]
            wij = weights[i, k]

            tot += x[i]*mu[j]*wij/2
            tot += x[j]*mu[i]*wij/2

    return tot

new_egs = []
for i in tqdm(range(counts.shape[0])):
    new_egs.append(
        new_eg(counts.values[i], mu, hs.neighbors.values, hs.weights.values)
    )

results['new_EG'] = new_egs
results['EG2'] = results.stdG**2 + results.EG**2
results['new_stdG'] = results['EG2']


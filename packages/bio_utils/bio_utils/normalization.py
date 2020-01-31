import numpy as np
import pandas as pd


def full_quantile(data):
    # Create a reference distribution

    ref = data.mean(axis=1).sort_values().reset_index()[0]
    ranks = data.rank(method='min').astype('int') - 1

    normed = {x: ranks[x].map(ref) for x in ranks.columns}
    normed = pd.DataFrame(normed)
    return normed

# DESEQ normalization/size factors


def deseq(data):
    row_geo_means = np.exp(np.log(data).mean(axis=1))
    ratios = data.divide(row_geo_means, axis=0)
    size_factors = ratios.median(axis=0)
    return size_factors


def nz_deseq(data):
    """
    Like the first, but only on the non-zero data
    """
    row_nz_count = (data > 1).sum(axis=1)
    row_geo_means = np.exp(np.log(data).sum(axis=1) / row_nz_count)
    ratios = data.divide(row_geo_means, axis=0)

    def map_fun(x):
        return x.loc[data[x.name] > 1].median()

    size_factors = ratios.apply(map_fun, axis=0)
    return size_factors

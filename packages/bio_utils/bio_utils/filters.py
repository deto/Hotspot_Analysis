# -*- coding: utf-8 -*-
"""Functions that are used to select genes

Functions here use a variety of criteria to reduce the number
of genes to a more manageable size - ideally extracting the
genes that are more biologically informative.

"""
from __future__ import absolute_import, print_function, division

import numpy as np
import pandas as pd


def filter_genes_threshold(data, threshold):
    """Filters out genes that are at least active in <threshold>
    number of samples.
    """

    num_detected = (data > 0).sum(axis=1)
    gene_passes = num_detected > threshold

    return gene_passes


def filter_genes_novar(data):
    """
    Filters out genes with 0 variance across samples
    :param data: Pandas DataFrame
    :return: Subset of data 0-variance rows removed
    """

    row_var = data.var(axis=1)
    gene_passes = row_var > 0

    return gene_passes


def filter_genes_fano(data, num_mad=2):
    """
    Uses fano filtering on the genes.  Splits into quantiles first.
    Only retains genes that have fano factor <num_mad> median absolute
        deviations in each quantile.

    :param data: pandas.dataframe
    :param num_mad: float
    :return: numpy.ndarray, bool, (Num_Genes)
        True for rows that pass the filter.  Otherwise, False;
    """
    genes = data.index
    data = data.values

    mu = np.mean(data, axis=1)
    sigma = np.std(data, axis=1)

    # slice up mu and sigma into bins, by mu
    aa = np.argsort(mu)
    mu_sort = mu[aa]
    sigma_sort = sigma[aa]

    N_QUANTS = 30
    m = mu_sort.size // N_QUANTS

    gene_passes = np.zeros(data.shape[0]) == 1

    for i in range(N_QUANTS):
        if(i == N_QUANTS - 1):
            rr = np.arange(i * m, mu_sort.size)
        else:
            rr = np.arange(i * m, (i + 1) * m)

        mu_quant = mu_sort[rr]
        mu_quant[mu_quant == 0] = 1  # so we don't divide by zero later
        sigma_quant = sigma_sort[rr]
        fano_quant = sigma_quant**2 / mu_quant
        mad_quant = np.median(np.abs(fano_quant - np.median(fano_quant)))
        gene_passes_quant = fano_quant > np.median(
            fano_quant) + num_mad * mad_quant
        gene_passes_quant_i = np.nonzero(gene_passes_quant)[0]
        gene_passes_i = gene_passes_quant_i + i * m
        gene_passes[gene_passes_i] = True

    # gene_passes are in sorted mu order, need to revert back to original order
    original_ii = np.argsort(aa)
    gene_passes = gene_passes[original_ii]

    return pd.Series(gene_passes, index=genes)

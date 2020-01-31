"""
Other enrichment methods
...requires other packages
"""
from __future__ import division, print_function
from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd
import tempfile
import os
import gseapy

from . import file_io
from .classes import GeneSet

def goseq_enrichment(selected_genes, bias_data, gene_sets, correct_false_neg=True):
    """
    Uses geoseq to calculated enrichment while accounting for the bias
    due to expression level

    selected_genes : list of str
        Genes which are DE

    bias_data : pandas.DataFrame
        Index : Gene ID
        Column1: Values to use for bias estimation

    gene_sets : GeneSet | list of GeneSet
        Sets to use for estimating enrichment

    correct_false_neg : bool
        Whether or not to use BH FDR procedure to correct p-vals in results
    """
    import goseq
    gs_file = tempfile.mktemp(prefix='goseq_')

    try:
        if(isinstance(gene_sets, GeneSet)):
            gene_sets = [gene_sets]

        file_io.write_gmt_file(gene_sets, gs_file)

        result = goseq.run(bias_data, selected_genes, gs_file)

    finally:
        os.remove(gs_file)

    # Change under_represented_pvalue to pvalue
    # Drop over_represented_pvalue

    result = result.drop("over_represented_pvalue", axis=1)

    cols = result.columns.tolist()
    ii = cols.index("under_represented_pvalue")
    cols[ii] = "pvalue"
    result.columns = cols

    # Correct False Negatives
    if correct_false_neg:

        reject, pv_corrected, aS, aB = multipletests(
            result['pvalue'], alpha=0.05, method='fdr_bh')
        result['pvalue'] = pv_corrected

    return result


def gsea(selected_genes, gene_sets, outdir):
    """
    Uses gseapy to calculate enrichment

    selected_genes : GeneSet
        Genes to test and values

    gene_sets : GeneSet | list of GeneSet
        Sets to use for estimating enrichment
    """

    assert isinstance(selected_genes, GeneSet)

    if selected_genes.values.std() == 0:
        raise ValueError("Selected genes either have no associated values or \
                         values are all identical")

    # Form the .rnk Data Frame
    rnk_df = pd.DataFrame({"Gene Name": selected_genes.values.index,
                           "Rank Stat": selected_genes.values.values})
    rnk_df = rnk_df.sort_values("Rank Stat", ascending=False)
    rnk_df = rnk_df.reset_index(drop=True)

    # Write gene_sets to a file (this is how gseapy expects them)
    gs_file = tempfile.mktemp(prefix='cara_', suffix=".gmt")

    try:
        if(isinstance(gene_sets, GeneSet)):
            gene_sets = [gene_sets]

        file_io.write_gmt_file(gene_sets, gs_file)

        result = gseapy.prerank(rnk_df, gs_file, outdir=outdir, max_size=1e6)

    finally:
        os.remove(gs_file)

    return result

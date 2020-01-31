"""Utility methods for computing statistics"""

from __future__ import division, print_function
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd
from tqdm import tqdm

from .classes import GeneSet

def all_v_all_enrichment(gene_sets_of_interest, gene_set_libraries, out_file, mode="overrepresentation", background=None, bias_data=None, progressbar=True):
    """
    Creates an Excel report with enrichment results
    Each library in gene_set_libraries gets a sheet
    In each sheet, there is a table showing enrichment results against each
    gene set in `gene_sets_of_interest`

    gene_sets_of_interest : list of GeneSet

    gene_set_libraries : dict of list of GeneSet
        key : Library name
        value: List of gene sets in the library

    out_file : str
        Name of file to create
        E.g. results.xlsx

    mode : str
        Either 'overrepresentation' or 'goseq'

    background : list of str
        Genes to use as a background when computing enrichments

    bias_data : pandas.DataFrame
        Index : Gene ID
        Column1: Values to use for bias estimation

    """

    if progressbar:
        total = len(gene_sets_of_interest) * len(gene_set_libraries)
        pbar = tqdm(total=total)

    if mode == "overrepresentation" and background is None:
        raise ValueError("Need background when mode is 'overrepresentation'")

    if mode == "goseq" and bias_data is None:
        raise ValueError("Need background when mode is 'goseq'")

    if mode == "goseq":
        from .other_analyses import goseq_enrichment

    enrichment_results = {}

    for es_name, es in gene_set_libraries.items():
        es_results = {}

        for gs in gene_sets_of_interest:
            selected_genes = gs.genes
            gene_sets = es

            if mode == "overrepresentation":
                result = gene_set_enrichment(
                    selected_genes, background, gene_sets)
            elif mode == "goseq":
                result = goseq_enrichment(selected_genes, bias_data, gene_sets)
            else:
                raise ValueError("Invalid option for argument: mode")

            result.index.name = gs.name

            es_results[gs.name] = result

            if progressbar:
                pbar.update()

        enrichment_results[es_name] = es_results

    writer = pd.ExcelWriter(out_file)
    for es_name, es in gene_set_libraries.items():
        sheetname = es_name
        offset = 0
        for gs in gene_sets_of_interest:
            result = enrichment_results[es_name][gs.name]
            result = result.loc[result.pvalue < 0.05]
            result.to_excel(writer, sheet_name=sheetname, startrow=offset)
            offset += result.shape[0] + 3

    writer.save()

    return enrichment_results


def gene_set_enrichment(selected_genes, all_genes, gene_sets, correct_false_neg=True):
    """
    Tests if each of the sets in gene_sets are enriched in the set
    of selected_genes, drawn from all_genes

    Genes in selected_genes and in each gene_set are tossed out prior to
    enrichment calculation if they are not in `all_genes`

    Parameters
    ----------
    selected_genes : Set of Strings or GeneSet
        Gene names to be tested for enrichment

    all_genes : List of Strings
        Gene names for all possible genes.  This is what
        selected_genes is drawn from

    gene_sets : GeneSet
        Sets of genes test agains the selected_genes

    correct_false_neg : bool, optional
        Whether or not to convert p-values to q-values using the BH procedure

    Returns
    -------
    result : pandas.DataFrame
        columns=["pvalue", "fc", "matches", "eff_set_size", "genes"])
        index=Gene Set Name for each gene set in `gene_sets`

    """
    if(isinstance(selected_genes, list)):
        selected_genes = set(selected_genes)

    if(isinstance(selected_genes, GeneSet)):
        selected_genes = selected_genes.genes

    # Cast to set for performance
    all_genes = set(all_genes)
    selected_genes = selected_genes & all_genes

    result = pd.DataFrame(index=[gs.name for gs in gene_sets],
                          columns=["pvalue", "fc", "matches", "eff_set_size", "genes"])

    for gs in gene_sets:

        gs_genes = gs.genes & all_genes
        if len(gs_genes) == 0:
            result.pvalue[gs.name] = 1.0
            result.fc[gs.name] = 1.0
            result.matches[gs.name] = 0
            result.eff_set_size[gs.name] = 0
            result.genes[gs.name] = ""
            continue

        selected_genes_in_set = len(gs_genes & selected_genes)

        result.pvalue[gs.name] = overrepresentation_unsigned(selected_genes,
                                                             all_genes, gs)

        actual_overlap = selected_genes_in_set / len(selected_genes)
        expected_overlap = len(gs_genes) / len(all_genes)

        result.fc[gs.name] = actual_overlap / expected_overlap
        result.matches[gs.name] = selected_genes_in_set
        result.eff_set_size[gs.name] = len(gs_genes)

        result.genes[gs.name] = ", ".join(sorted(gs_genes & selected_genes))

    if(correct_false_neg):
        reject, pv_corrected, aS, aB = multipletests(
            result['pvalue'].values, alpha=0.05, method='fdr_bh')

        result['pvalue'] = pv_corrected

    result.sort_values("pvalue", inplace=True)

    return result


def overrepresentation_unsigned(selected_genes, all_genes, gene_set):
    """
    Compute whether `selected_genes` is overrepresented in `gene_set` given
    the background set of `all_genes`
    """
    gs_genes = gene_set.genes
    selected_genes_in_set = len(gs_genes & selected_genes)
    gs_genes_in_all = len(gs_genes & all_genes)

    M = len(all_genes)
    N = len(selected_genes)

    p_val = hypergeom.sf(selected_genes_in_set - 1, M,
                         gs_genes_in_all, N)

    return p_val


def gs_jaccard(gs_a, gs_b):
    """
    Computes the jaccard similarity between two gene sets

    Parameters
    ----------
    gs_a : GeneSet
        First gene set
    gs_b : GeneSet
        Second gene set

    Returns
    -------
    float
        Jaccard similarity index (from 0 to 1.0)
    """

    s1 = gs_a.genes
    s2 = gs_b.genes

    s1_intersect_s2 = len(s1 & s2)

    if(s1_intersect_s2 == 0):
        return 0.0

    return s1_intersect_s2 / (len(s1) + len(s2) - s1_intersect_s2)

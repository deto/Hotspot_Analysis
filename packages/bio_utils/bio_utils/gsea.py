import numpy as np
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
from gseapy.algorithm import enrichment_score, gsea_significance


def gsea_signed(values, signatures, midpoint=0, processes=4,
                weighted_score_type=0, progressbar=True):
    """
    Compute the signed GSEA by flipping the rankings for negative genes
    prior to running regular GSEA

    values: pandas.Series
        values to use for ranking and index is gene names

    signatures: list of gene_enrich.GeneSet

    """

    # For each gs in GMT

    enrichment_scores = []
    enrichment_nulls = []
    hit_ind = []
    rank_ES = []
    LE_genes = []
    names = [x.name for x in signatures]

    pool = Pool(processes=processes)

    partial_map_fun = partial(
        gsea_signed_single, values=values, midpoint=midpoint,
        weighted_score_type=weighted_score_type
    )

    if progressbar:
        pbar = tqdm(total=len(signatures))

    for (es, ind, RES, bg, le_g) in pool.imap(partial_map_fun, signatures):
        # es, ind, RES, bg = gsea_signed_single(values, sig, midpoint)
        enrichment_scores.append(es)
        rank_ES.append(RES)
        hit_ind.append(ind)
        enrichment_nulls.append(bg)
        LE_genes.append(le_g)

        if progressbar:
            pbar.update()

    if progressbar:
        pbar.close()

    pool.close()
    pool.join()

    sig_table = gsea_significance(enrichment_scores, enrichment_nulls)
    sig_table_df = pd.DataFrame(list(sig_table), index=names,
                                columns=['ES', 'nES', 'pval', 'fdr'])

    gsea_details = {names[i]: (hit_ind[i], rank_ES[i], LE_genes[i])
                    for i in range(len(names))}

    # Add leading-edge genes to the table
    sig_table_df['Leading-Edge Genes'] = \
        [", ".join(gsea_details[i][2]) for i in sig_table_df.index]

    return sig_table_df, gsea_details


def gsea_signed_single(sig, values, midpoint=0, weighted_score_type=0):
    """
    Compute signed GSEA on a single signature
    """

    assert values.index.is_unique, "Index for values must be unique"

    v2 = values.copy()
    for gene in sig.values.index:
        if sig.values[gene] < 0 and gene in v2.index:
            v2[gene] = (v2[gene]-midpoint)*-1 + midpoint

    # convert signature from signed into gseapy format
    v2 = v2.sort_values(ascending=False)
    gene_list = v2.index.tolist()
    correl_vector = v2.values
    sig_genes = list(sig.genes)
    rs = np.random.RandomState(1028)

    # gsea
    es, esnull, ind, RES = enrichment_score(
        gene_list, gene_set=sig_genes,
        correl_vector=correl_vector,
        weighted_score_type=weighted_score_type,
        nperm=1000, rs=rs)

    # Add in the leading-edge genes
    LE_genes = []
    if es > 0:
        ii = np.argmax(RES)
        for hi in ind:
            if hi > ii:
                break
            gene = gene_list[hi]
            if sig.values[gene] < 0:
                gene = "(" + gene + ")"
            LE_genes.append(gene)
    else:
        ii = np.argmin(RES)
        for hi in ind[::-1]:
            if hi < ii:
                break
            gene = gene_list[hi]
            if sig.values[gene] < 0:
                gene = "(" + gene + ")"
            LE_genes.append(gene)

    return es, ind, RES, esnull, LE_genes

from __future__ import division, print_function
import json
import pandas as pd

from .orthology import get_human_mouse_ortho
from .classes import GeneSet


def load_gene_set_gmt(file_name, mouse2human=False, human2mouse=False):
    """
    Loads the gene sets in the gmt file <file_name>

    Parameters
    ==========
    file_name : str
        Name of file with the gene set
    mouse2human : bool
        If True, use the MGI orthologs to convert mouse genes to human orthologs
    human2mouse : bool
        If True, use the MGI orthologs to convert human genes to mouse orthologs

    returns: List of GeneSet's
    """

    plus_terms = {'plus', 'up'}
    minus_terms = {'minus', 'dn', 'down'}

    out_gene_sets = {}
    with open(file_name, 'r') as fin:
        for line in fin:
            line_s = line.strip().split('\t')
            if len(line_s) <= 2:
                continue  # Skip bad lines
            set_name = line_s[0]
            genevals = line_s[2:]  # Skip the comment portion

            # Remove trailing _up or _dn from set_name
            set_sign = 1.0
            name_split = set_name.split("_")

            if name_split[-1].lower() in plus_terms:
                set_name = '_'.join(name_split[:-1])

            elif name_split[-1].lower() in minus_terms:
                set_name = '_'.join(name_split[:-1])
                set_sign = -1.0

            # Genes can be stored with values
            # E.g. HNF4,1.5
            # Here, split on the comma
            genes = []
            values = []
            for gv in genevals:
                gv_s = gv.split(",")
                if len(gv_s) == 1:
                    genes.append(gv_s[0])
                    values.append(set_sign)
                elif len(gv_s) == 2:
                    genes.append(gv_s[0])
                    values.append(float(gv_s[1]))
                else:
                    raise Exception("Error on line: " + line[0:20])

            # Normalize all values to upper-case
            genes = [x.upper() for x in genes]

            description = line_s[1]

            if set_name not in out_gene_sets:
                new_set = GeneSet(set_name, genes, description,
                                  file_name, values=values)
                out_gene_sets[set_name] = new_set
            else:
                # Just update the existing set with new values
                existing_set = out_gene_sets[set_name]
                existing_set.add_genes(genes, values)

    out_gene_sets = list(out_gene_sets.values())

    if mouse2human:
        _ortho2human, _ortho2mouse, _human2ortho, _mouse2ortho = get_human_mouse_ortho()
        for gs in out_gene_sets:
            hgenes = []
            hvals = []
            for mgene, mval in gs.values.iteritems():
                for ortho_id in _mouse2ortho[mgene]:
                    for hg in _ortho2human[ortho_id]:
                        hgenes.append(hg)
                        hvals.append(mval)

            gs.genes = set(hgenes)
            gs.values = pd.Series(
                {x: y for x, y in zip(hgenes, hvals)}
            )

    if human2mouse:
        _ortho2human, _ortho2mouse, _human2ortho, _mouse2ortho = get_human_mouse_ortho()
        for gs in out_gene_sets:
            mgenes = []
            mvals = []
            for hgene, hval in gs.values.iteritems():
                for ortho_id in _human2ortho[hgene]:
                    for mg in _ortho2mouse[ortho_id]:
                        mgenes.append(mg)
                        mvals.append(hval)

            gs.genes = set(mgenes)
            gs.values = pd.Series(
                {x: y for x, y in zip(mgenes, mvals)}
            )

    return out_gene_sets


def write_gmt_file(gene_sets, file_name, include_values=False):
    """Writes a .gmt file with the sets in `gene_sets`

    Parameters
    ----------
    gene_sets : list of GeneSet
        Sets to write to file
    file_name : str
        File to write

    """
    with open(file_name, 'w') as fout:
        for gs in gene_sets:

            if include_values:
                gene_vals = [g + "," + str(v)
                             for g, v in gs.values.iteritems()]
                fout.write(
                    "\t".join([gs.name, gs.description] + gene_vals) + "\n")
            else:
                fout.write(
                    "\t".join([gs.name, gs.description] + list(gs.genes)) + "\n")


def load_gene_set_json(filename):
    """
    Load a list of gene sets from the structured
    json signature format
    """
    with open(filename, 'r') as fin:
        data = json.load(fin)

    sig_set = []

    for sig in data:
        gs = GeneSet(sig['name'], sig['values'].keys(),
                     sig['description'], file_name=filename,
                     values=sig['values'].values())
        sig_set.append(gs)

    return sig_set

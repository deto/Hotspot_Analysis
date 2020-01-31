"""
A Module for working with gene sets.
- Usually stored in .gmt files
"""

from __future__ import division, print_function
import pandas as pd


class GeneSet:
    """
    Wrapper class for sets of genes
    """

    def __init__(self, name, genes, description="", file_name='None', values=None):
        """
        Set of Genes in the set (Set of String)

        name: str
        genes: list of str
        description: str
        file_name: str
        values: list of float
        """

        self.name = name

        self.description = description

        self.genes = {x.upper() for x in genes}

        self.file_name = file_name  # Source File for set

        # Values hold the 'sign' information on genes
        if values is None:
            self.values = pd.Series(1.0, index=self.genes)
            self.signed = False
        else:
            assert len(values) == len(genes)
            g_dict = {x.upper(): y for x,y in zip(genes, values)}
            self.values = pd.Series(g_dict)
            self.signed = True

    def add_genes(self, genes, values=None):
        if values is None:
            values = [1.0 for x in genes]

        genes = [x.upper() for x in genes]
        g_dict = {x: y for x, y in zip(genes, values)}
        new_values = pd.Series(g_dict)

        self.values = self.values.append(new_values)

        to_drop = self.values.index.duplicated(keep='last')
        self.values = self.values[~to_drop]

        self.genes |= set(genes)

    def __str__(self):
        return "GeneSet: '" + self.name + "' with " + str(len(self.genes)) + " Genes"

    def __repr__(self):
        return "<" + str(self) + ">"

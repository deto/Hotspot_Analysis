# gene_enrich
Enrichment Analysis on Gene Sets

# Installation
1. Clone repository
2. In the repository's top-level directory (with setup.py) run:

    pip install -e .
    
When installed this way, it is easy to update the library if changes have been pushed to Github.

To update, just run `git pull` in the repository.

# Load gene sets from .gmt files
```python
import gene_enrich
gene_sets = gene_enrich.load_gene_set_gmt("/path/to/your/gmt/file/H_Hallmark.gmt")
print(gs[0])
# GeneSet: 'HALLMARK_ADIPOGENESIS' with 200 Genes
```

`gene_sets` is a list of **GeneSet** - a class which encapsultes the name, genes and values in a gene set.

```python
print(gs[0].name)
# HALLMARK_ADIPOGENESIS

print(gs[0].genes)
# set(['RNF11', 'NDUFB7', 'CD36', 'DHCR7', ... , 'SLC1A5', 'SDHC', 'BCL6', 'IMMT', 'COQ3'])

print(gs[0].values)
"""
ABCA1       1.0
ABCB8       1.0
ACAA2       1.0
ACADL       1.0
ACADM       1.0
ACADS       1.0
ACLY        1.0
"""
```

# Where to download gene sets (.gmt files)
1. [MSIGDB - http://software.broadinstitute.org/gsea/downloads.jsp](http://software.broadinstitute.org/gsea/downloads.jsp)
2. [Enrichr - http://amp.pharm.mssm.edu/Enrichr/#stats](http://amp.pharm.mssm.edu/Enrichr/#stats)

# Running an enrichment analysis (one vs many)
For example, if you have a list of differential genes and want to see if they are enriched for any Hallmark gene sets

```python
import gene_enrich

# Load the Hallmark set
hallmark_sets = gene_enrich.load_gene_set_gmt("/path/to/your/gmt/file/H_Hallmark.gmt")

# Put your differential genes into a set
diff_genes = set(['PCK', 'GAPDH', 'GCK', 'SLC2A2' ... ])

# Define a background - this should be all genes that could have been called as differential
background = set(['A2M, 'AAAS', 'AADAT', 'AARS', 'ABAT', 'ABCA1', ..., 'ZPBP', 'ZW10', 'ZWINT', 'ZYX'])

# Run the analysis
result = gene_enrich.gene_set_enrichment(diff_genes, background, hallmark_sets)

"""Result:
                                         pvalue       fc matches eff_set_size
HALLMARK_COMPLEMENT                           0    21.93     200          200   
HALLMARK_COAGULATION                3.69318e-41  9.05804      57          138   
HALLMARK_INTERFERON_GAMMA_RESPONSE  0.000486657  2.52195      23          200   
HALLMARK_INFLAMMATORY_RESPONSE       0.00308205  2.30265      21          200   
HALLMARK_APOPTOSIS                     0.142707  1.90696      14          161   

                                                                                genes  
HALLMARK_COMPLEMENT                 ACTN2, ADAM9, ADRA2B, AKAP10, ANG, ANXA5, APOA...  
HALLMARK_COAGULATION                ADAM9, ANG, APOC1, C1QA, C1R, C1S, C2, C3, C9,...  
HALLMARK_INTERFERON_GAMMA_RESPONSE  C1R, C1S, CASP1, CASP3, CASP4, CASP7, CCL5, CF...  
HALLMARK_INFLAMMATORY_RESPONSE      CCL5, CD55, F3, GNAI3, GP1BA, IL6, IRF1, IRF7,...  
HALLMARK_APOPTOSIS                  CASP1, CASP3, CASP4, CASP7, CASP9, CLU, F2, IL...  
"""
```

# Running an analysis of many gene sets vs many libraries (many vs many)
For running an analysis of many sets of genes against many gene set libraries
Results are output in a formatted excel file and a dictionary data structure

```python
# Load gene set libraries into a dictionary
gs_libraries = {
  'hallmark': gene_enrich.load_gene_set_gmt("/path/to/your/gmt/file/H_Hallmark.gmt"), 
   'GO': gene_enrich.load_gene_set_gmt("/path/to/your/gmt/file/C5_GO_ALL.gmt"), 
  'KEGG': gene_enrich.load_gene_set_gmt("/path/to/your/gmt/file/C5_Kegg.gmt")
}

# Wrap your gene sets in GeneSet classes
diff_genes1 = gene_enrich.GeneSet('Set #1', genes=['FASN', 'CYP3A4', 'CD36'...])
diff_genes2 = gene_enrich.GeneSet('Set #1', genes=['SIRT2', 'MAPK', 'CYP2B6'...])
diff_genes3 = gene_enrich.GeneSet('Set #1', genes=['IL6', 'CCL5', 'CASP3'...])
my_sets = [diff_genes1, diff_genes2, diff_genes3]
 
# Define a background as before
background = set(['A2M, 'AAAS', 'AADAT', 'AARS', 'ABAT', 'ABCA1', ..., 'ZPBP', 'ZW10', 'ZWINT', 'ZYX'])

# Run the analyses
result = gene_enrich.all_v_all_enrichment(my_sets, gs_libraries, out_file="output.xlsx", background=background)

# Note: background must be provided as a keyword argument
# Outputs are written to out_file.  Only enrichments with p<0.05 are kept
# result contains a hierarchical dictionary with all the results
```

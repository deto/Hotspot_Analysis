import ez_setup
ez_setup.use_setuptools()


from setuptools import setup, find_packages

setup(
    name="gene_enrich",
    version="0.0.1",
    packages=find_packages(),

    package_data={
        '': ['MGI/*']
    },

    install_requires=["pandas>=0.19.0", "statsmodels>=0.6.1",
                      "tqdm>=4.8.2", "scipy>=0.18.0"],

    author="David DeTomaso",
    author_email="David.DeTomaso@berkeley.edu",
    description="Package for enrichment analysis on gene sets",
    keywords="Bioinformatics genomics rna-seq",
)

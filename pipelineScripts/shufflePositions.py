"""
Shuffles the positions for barcodes -> used to compare null model assumptions
"""
import pandas as pd

positions_in = snakemake.input['positions']
positions_out = snakemake.output['positions']

positions = pd.read_table(positions_in, index_col=0)

positions_shuff = positions.sample(frac=1)
positions_shuff.index = positions.index

positions_shuff.to_csv(positions_out, sep='\t')

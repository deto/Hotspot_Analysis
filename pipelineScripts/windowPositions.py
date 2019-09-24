"""
Shuffles the positions for barcodes -> used to compare null model assumptions
"""
import pandas as pd

positions_in = snakemake.input['positions']
positions_out = snakemake.output['positions']

xmin = snakemake.params['xmin']
xmax = snakemake.params['xmax']
ymin = snakemake.params['ymin']
ymax = snakemake.params['ymax']

positions = pd.read_table(positions_in, index_col=0)

positions = positions.loc[
    (positions.iloc[:, 0] >= xmin) &
    (positions.iloc[:, 0] <= xmax) &
    (positions.iloc[:, 1] >= ymin) &
    (positions.iloc[:, 1] <= ymax)
]

positions.to_csv(positions_out, sep='\t')

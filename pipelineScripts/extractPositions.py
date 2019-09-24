"""
This script just extracts the Positions in the loompy file

Mainly used so that it can be substituted as an endpoint in pipelines
where the positions may also be used as the latent space
"""
import pandas as pd
import loompy

loom_file = snakemake.input['loom']
positions_out = snakemake.output['positions']

with loompy.connect(loom_file, 'r') as ds:
    positions = ds.ca['Position'][:]
    barcodes = ds.ca['Barcode'][:]


positions = pd.DataFrame(
    positions,
    index=barcodes,
    columns=['Comp{}'.format(i+1) for i in range(positions.shape[1])]
)
positions.to_csv(positions_out, sep='\t')

import pandas as pd
from sklearn.manifold import TSNE
import sys

latent_file = snakemake.input['latent']
out_file = snakemake.output['out']

pc_data = pd.read_table(latent_file, index_col=0)


model = TSNE(n_components=2, verbose=1)

result = model.fit_transform(pc_data)

result = pd.DataFrame(result, index=pc_data.index, columns=['tsne1', 'tsne2'])

result.to_csv(out_file, sep='\t')

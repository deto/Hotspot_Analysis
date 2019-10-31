import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loompy


hs_results = pd.read_table(
    "../hotspot/hotspot_shuffled.txt", index_col=0
)

ds = loompy.connect("../../data/10x_PBMC_w_proteins/cd4/data.loom", "r")
gene_info = pd.DataFrame({
    "EnsID": ds.ra["EnsID"][:],
    "Symbol": ds.ra["Symbol"][:]
}).set_index("EnsID")
ens_map = {k: v for k, v in gene_info.itertuples()}
counts = pd.DataFrame(ds.layers['scaled'][:, :], index=gene_info.index)
ds.close()

gene_mean = counts.mean(axis=1)


# %%
# Distribution looks pretty good actually

fig, axs = plt.subplots(1, 2, figsize=(9, 4))

plt.sca(axs[0])
plt.hist(hs_results.Z,  50)
plt.xlabel('Z-score')

plt.sca(axs[1])
data = hs_results.join(gene_mean.rename("Mean"))
plt.plot(data.Mean, data.Z, 'o', ms=2)
plt.xscale('log')
plt.xlabel('Gene Mean')
plt.ylabel('Z-score')

plt.suptitle('Hotspot - CD4 Shuffled')
plt.show()

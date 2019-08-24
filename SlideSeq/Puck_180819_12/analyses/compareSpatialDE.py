import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bio_utils import hover_plot



sp_res = pd.read_table("../spatialDE/spatialDE.txt", index_col=0)
sp_res.head()

# Observation 1 - all the values are significant
#   Over 98 percent significant at FDR < 0.05


# What if we do the correction? They use a strange-looking custom script
from statsmodels.stats.multitest import multipletests
q2 = multipletests(sp_res.pval, method='fdr_bh')[1]
# It's a little different, but still get 5789 candidates


hs_res = pd.read_table("../hotspot/hotspot.txt", index_col=0)

# Only plot common variables

common = sp_res.index & hs_res.index

plt.figure()

plt.plot(
    sp_res.LLR[common],
    hs_res.Z[common],
    'o', ms=2
)
plt.show()



plt.figure()

plt.plot(
    np.log10(sp_res.pval[common]),
    np.log10(hs_res.Pval[common]),
    'o', ms=2
)
plt.show()



hover_plot(
    np.log10(sp_res.pval[common]+1e-100),
    np.log10(hs_res.Pval[common]+1e-100),
    common,
    'o', ms=2
)
plt.show()

# Maybe plot ranking vs ranking?

p1 = sp_res.pval[common].sort_values()
r1 = pd.Series(np.arange(len(common))+1, index=p1.index)

p2 = hs_res.Pval[common].sort_values()
r2 = pd.Series(np.arange(len(common))+1, index=p2.index)


data = pd.concat(
    (r1.rename('sDE-rank'), r2.rename('hs-rank')), axis=1)

data = data.loc[data.min(axis=1) <= 600]
data = data.join(sp_res.l)

hover_plot(
    data['sDE-rank'],
    data['hs-rank'],
    data.index,
    'o', ms=2
)
plt.xlabel('sDE-rank')
plt.ylabel('hs-rank')
plt.show()


# What overlap is their in the top N?
N = 1000
s1 = sp_res.pval[common].sort_values().index[0:N]
s2 = hs_res.Pval[common].sort_values().index[0:N]

print(len(s1 & s2))
print(len(s1 | s2))

# Does the length-scale have something to do with it?
# HS will miss really long-range correlations as we limited to K=30 for this trial

plt.figure()

plt.scatter(
    x=data['sDE-rank'],
    y=data['hs-rank'],
    c=np.log10(data['l']),
    s=2, cmap='jet'
)
plt.colorbar()
plt.show()

# Kind of
# - insignificant in sDE but sig in HS is usually the lowest length scale
# - however, at .25 this length is smaller than the distance between barcodes
#   and almost meaningless here

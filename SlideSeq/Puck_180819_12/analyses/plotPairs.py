import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


Z = pd.read_table("../hotspot/hotspot_pairs_z.txt.gz", index_col=0)
hs_results = pd.read_table("../hotspot/hotspot.txt", index_col=0)

# Cluster things!

# to_drop = ['mt-Rnr1', 'mt-Rnr2']
# Z = Z.drop(to_drop, axis=1).drop(to_drop, axis=0)

sns.clustermap(Z, vmin=-2, vmax=2,
               metric='correlation', yticklabels=True, method='average')
plt.show()

# Now cluster and divide

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram


def sort_clusters(cl):
    map_fun = {old_i: new_i+1 for new_i, old_i in enumerate(cl.value_counts().index)}
    cl = cl.map(map_fun)
    return cl


L = linkage(Z, method='average', metric='correlation')

cl = fcluster(L, t=20, criterion='maxclust')
cl = pd.Series(cl, index=Z.index)
cl = sort_clusters(cl)


cmap = plt.get_cmap('tab20')

cluster_colors = [cmap.colors[i % len(cmap.colors)] for i in cl]

sns.clustermap(Z, vmin=-2, vmax=2,
               metric='correlation', yticklabels=True, row_colors=cluster_colors)
plt.show()


m1 = cl.index[cl == 1]
m2 = cl.index[cl == 2]
m3 = cl.index[cl == 3]
m4 = cl.index[cl == 4]
m5 = cl.index[cl == 5]

Z_sub = Z.loc[m1, m1]
Z_sub = Z_sub.loc[hs_results.loc[m1].Z.sort_values().index]
sns.clustermap(Z_sub, vmin=-2, vmax=2,
               metric='correlation', yticklabels=True, method='average',
               row_cluster=False)
plt.show()

L = linkage(Z_sub, method='average', metric='correlation')

cl_sub = fcluster(L, t=5, criterion='maxclust')
cl_sub = pd.Series(cl_sub, index=Z_sub.index)
cl_sub = sort_clusters(cl_sub)

cmap = plt.get_cmap('tab20')

cluster_colors = [cmap.colors[i % len(cmap.colors)] for i in cl_sub]

sns.clustermap(Z_sub, vmin=-2, vmax=2,
               metric='correlation', yticklabels=True, row_colors=cluster_colors)
plt.show()


hs_results.loc[m1].sort_values('Z')


# How do the spatialDE results look here?
sp_res = pd.read_table("../spatialDE/spatialDE.txt", index_col=0)


selected = sp_res.sort_values('LLR').index[-500:]

sp_colors = ['black' if x in selected else 'white' for x in Z.index]
sns.clustermap(Z, vmin=-2, vmax=2,
               metric='correlation', yticklabels=True, row_colors=sp_colors)
plt.show()



sp_colors = ['black' if x in selected else 'white' for x in Z_sub.index]
sns.clustermap(Z_sub, vmin=-2, vmax=2,
               metric='correlation', yticklabels=True, row_colors=sp_colors)
plt.show()


in_sp = pd.Series([1 if x in selected else 0 for x in hs_results.index],
                  index=hs_results.index, name='InSp')

hs_results.join(in_sp).loc[m4].sort_values('Z').tail(40)

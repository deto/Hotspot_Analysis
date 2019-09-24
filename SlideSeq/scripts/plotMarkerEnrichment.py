# %% Setup

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve

# %% Load our different orderings

scores = {}

hs_file = "../Puck_180819_12/hotspot/hotspot_300.txt"
hs_results = pd.read_table(hs_file, index_col=0)
scores["Hotspot-Position"] = hs_results.Z

hs_file = "../Puck_180819_12/hotspot/hotspot_pca_threshold.txt"
hs_results = pd.read_table(hs_file, index_col=0)
scores["Hotspot-PCA-Thresh"] = hs_results.Z

hs_file = "../Puck_180819_12/hotspot/hotspot_pca_hvg.txt"
hs_results = pd.read_table(hs_file, index_col=0)
scores["Hotspot-PCA-Hvg"] = hs_results.Z

cluster_file = "../Puck_180819_12/clusters/markers_threshold.txt.gz"
clusters = pd.read_table(cluster_file)
cl_p = clusters.pivot(index="gene", columns="cluster", values="LR")
scores["DE-PCA-Thresh"] = cl_p.max(axis=1)

cluster_file = "../Puck_180819_12/clusters/markers_hvg.txt.gz"
clusters = pd.read_table(cluster_file)
cl_p = clusters.pivot(index="gene", columns="cluster", values="LR")
scores["DE-PCA-Hvg"] = cl_p.max(axis=1)


# %% Load Markers

markers_file = "../DropViz/edger_markers_1vAll.txt"
pos_data = pd.read_table(markers_file)

N = 100  # Number of genes per cluster

markers = {}
for c_id, group in pos_data.groupby("Cluster"):
    group = group.loc[group.FDR < 0.1]
    m = group.sort_values("LR", ascending=False).GeneSymbol[0:N]
    m = set(m)
    markers[c_id] = m

# Collapse markers
from functools import reduce

markers_all = reduce(set.union, markers.values(), set())
markers_all = {x.lower() for x in markers_all}

# Load manual markers
markers_file = "../DropViz/manual_markers.txt"
pos_data = pd.read_table(markers_file, comment="#", header=None)
markers_manual = {x.lower() for x in pos_data.iloc[:, 0]}

# %% Plot ROC


def calc_roc(scores, markers):
    gene_true = [x.lower() in markers for x in scores.index]
    fpr, tpr, thresholds = roc_curve(gene_true, scores.values)
    return fpr, tpr


def calc_pr(scores, markers):
    gene_true = [1 if x.lower() in markers else 0 for x in scores.index]
    precision, recall, thresholds = precision_recall_curve(
        gene_true, scores.values
    )
    return precision, recall


fig, axs = plt.subplots(2, 2, figsize=(8, 8))


methods = list(scores.keys())
markers = markers_all

plt.sca(axs[0, 0])
for method in methods:
    fpr, tpr = calc_roc(scores[method], markers)
    plt.plot(fpr, tpr, label=method)

plt.title('ROC - Computational Markers')

plt.sca(axs[0, 1])
for method in methods:
    precision, recall = calc_pr(scores[method], markers)
    plt.plot(recall, precision, label=method)

plt.title('PR - Computational Markers')

markers = markers_manual

plt.sca(axs[1, 0])
for method in methods:
    fpr, tpr = calc_roc(scores[method], markers)
    plt.plot(fpr, tpr, label=method)

plt.title('ROC - Manual Markers')

plt.sca(axs[1, 1])
for method in methods:
    precision, recall = calc_pr(scores[method], markers)
    plt.plot(recall, precision, label=method)

plt.title('PR - Manual Markers')

plt.legend()
plt.subplots_adjust(hspace=.4)
plt.show()
# plt.savefig('ROC_PR.png', dpi=300)

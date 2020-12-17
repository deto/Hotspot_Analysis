from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc
from functools import reduce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import loompy

plt.rcParams['svg.fonttype'] = 'none'

# %% First get the markers from DropViz data

marker_de = pd.read_table("../../DropViz/ranksums_markers_1vAll.txt")

groups = {
    l: marker_de.loc[marker_de.level == l].set_index('gene')
    for l in marker_de.level.unique()
}


def getMarkers(df, N):

    markers = df.loc[df.FDR < .1].sort_values('AUC').index.tolist()
    if len(markers) > N:
        markers = markers[:N]

    return markers


N = 100
markers = {l: getMarkers(v, N) for l, v in groups.items()}

all_markers = reduce(
    lambda x, y: x | set(y), markers.values(), set()
)

all_markers_norm = set([x.upper() for x in all_markers])

print(len(all_markers))
print(len(all_markers_norm))


# %% Define some lists...

k_values = [5, 10, 30, 50, 100, 300, 500, 1000]

hs_files = {
    k: "../../Puck_180819_12/k_sensitivity/k_{k}/hotspot.txt".format(k=k) for k in k_values
}


# %% Parse data for AUC and ROC curves

aucs = []

for k in hs_files.keys():
    ff = hs_files[k]
    res = pd.read_table(ff, index_col=0)

    y_score = res.Z.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    auc_i = roc_auc_score(y_true, y_score)

    aucs.append([k, auc_i])

aucs = pd.DataFrame(aucs, columns=['K', 'AUROC'])

# %%

plt.figure()
ax = plt.gca()
plt.plot(aucs['K'], aucs['AUROC'], 'o-')
plt.xscale('log')
plt.xticks(aucs['K'].unique())
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_ticks([], minor=True)
plt.xlabel('# of Neighbors')
plt.ylabel('Area Under ROC Curve')
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', ls=(0, (5, 5)), lw=.5)
plt.savefig('AU_ROC_k_sensitivity.svg', dpi=300)

# %% What about the precision recall curve?

auprs = []

for k in hs_files.keys():
    ff = hs_files[k]
    res = pd.read_table(ff, index_col=0)

    y_score = res.Z.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    aupr = auc(recall, precision)
    auprs.append([k, aupr])


auprs = pd.DataFrame(auprs, columns=['K', 'AUPRC'])

# %%

plt.figure()
ax = plt.gca()
plt.plot(auprs['K'], auprs['AUPRC'], 'o-')
plt.xscale('log')
plt.xticks(auprs['K'].unique())
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_ticks([], minor=True)
plt.ylabel('Area Under\nPrecision-Recall Curve')
plt.xlabel('# of Neighbors')
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', ls=(0, (5, 5)), lw=.5)
plt.savefig('AU_PR_k_sensitivity.svg', dpi=300)

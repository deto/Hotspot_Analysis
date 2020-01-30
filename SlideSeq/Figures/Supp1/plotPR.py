from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc
from functools import reduce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import loompy


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


hs_files = [
    "../../Puck_180819_9/hotspot/hotspot.txt",
    "../../Puck_180819_10/hotspot/hotspot.txt",
    "../../Puck_180819_11/hotspot/hotspot.txt",
    "../../Puck_180819_12/hotspot/hotspot.txt",
]

sde_files = [
    "../../Puck_180819_9/spatialDE/spatialDE.txt",
    "../../Puck_180819_10/spatialDE/spatialDE.txt",
    "../../Puck_180819_11/spatialDE/spatialDE.txt",
    "../../Puck_180819_12/spatialDE/spatialDE.txt",
]

# %% Make a ROC Curve

aucs = []

plt.figure()
for ff in hs_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.Z.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_i = roc_auc_score(y_true, y_score)
    aucs.append([auc_i, 'Hotspot'])

    plt.plot(fpr, tpr, color='blue')

for ff in sde_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.LLR.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_i = roc_auc_score(y_true, y_score)
    aucs.append([auc_i, 'spatialDE'])

    plt.plot(fpr, tpr, color='green')

plt.show()


aucs = pd.DataFrame(aucs, columns=['AUC', 'Method'])

sns.boxplot(x='Method', y='AUC', data=aucs)
plt.ylim(0, 1)
plt.show()

# %% Ok - this is very promising.  What if we just looked at average expression tho...
pucks = [
    "Puck_180819_9",
    "Puck_180819_10",
    "Puck_180819_11",
    "Puck_180819_12",
]

puck_means = {}
for puck in pucks:
    ls = loompy.connect("../../../data/SlideSeq/{}/data.loom".format(puck), mode="r")
    counts = ls.layers['scaled'][:, :]
    gene_info = ls.ra['EnsID', 'Symbol']
    ls.close()

    # Have to do this because data_slideseq makes it a numpy array
    gene_info = pd.DataFrame(
        gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
    counts = pd.DataFrame(counts, index=gene_info.index)

    # need counts, latent, and num_umi
    valid_genes = (counts > 0).sum(axis=1) > 50
    gene_means = counts.loc[valid_genes].mean(axis=1)
    puck_means[puck] = gene_means
    del counts

# %%


aucs = []

plt.figure()

for ff in hs_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.Z.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_i = roc_auc_score(y_true, y_score)
    aucs.append([auc_i, 'Hotspot'])

    if ff == hs_files[-1]:
        plt.plot(fpr, tpr, color='blue', label="Hotspot")
    else:
        plt.plot(fpr, tpr, color='blue')

for ff in sde_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.LLR.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_i = roc_auc_score(y_true, y_score)
    aucs.append([auc_i, 'spatialDE'])

    if ff == sde_files[-1]:
        plt.plot(fpr, tpr, color='green', label="spatialDE")
    else:
        plt.plot(fpr, tpr, color='green')

i = 0
for puck in puck_means:
    y_score = puck_means[puck].values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in puck_means[puck].index]
    ).astype('int')

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_i = roc_auc_score(y_true, y_score)
    aucs.append([auc_i, 'Expression'])

    if i == 0:
        plt.plot(fpr, tpr, color='orange', label='Mean Expression')
        i = 1
    else:
        plt.plot(fpr, tpr, color='orange')


plt.legend()
plt.show()

aucs = pd.DataFrame(aucs, columns=['AUC', 'Method'])

sns.stripplot(x='Method', y='AUC', data=aucs)
plt.ylim(0, 1)
plt.show()


# %% What about the precision recall curve?

auprs = []

plt.figure()

for ff in hs_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.Z.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    aupr = auc(recall, precision)
    auprs.append([aupr, 'Hotspot'])

    plt.plot(recall, precision, color='blue')

for ff in sde_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.LLR.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    aupr = auc(recall, precision)
    auprs.append([aupr, 'spatialDE'])

    plt.plot(recall, precision, color='green')

for puck in puck_means:
    y_score = puck_means[puck].values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in puck_means[puck].index]
    ).astype('int')

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    aupr = auc(recall, precision)
    auprs.append([aupr, 'Expression'])

    plt.plot(recall, precision, color='orange')

plt.show()

auprs = pd.DataFrame(auprs, columns=['AUPR', 'Method'])

sns.stripplot(x='Method', y='AUPR', data=auprs)
#plt.ylim(0, 1)
plt.show()

# Mean-method has much lower precision
# Not going to use to final figure, but good to know

# %% Now create the final figure

hs_color = sns.color_palette("deep")[0]
sde_color = sns.color_palette("deep")[2]

auprs = []

fig, axs = plt.subplots(
    1, 2,
    gridspec_kw={
        'width_ratios': [2, 1],
        'wspace': .5,
    }
)
plt.sca(axs[0])

for ff in hs_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.Z.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    aupr = auc(recall, precision)
    auprs.append([aupr, 'Hotspot'])

    label = 'Hotspot' if ff == hs_files[-1] else None
    plt.plot(recall, precision, color=hs_color, label=label)

for ff in sde_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.LLR.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    aupr = auc(recall, precision)
    auprs.append([aupr, 'spatialDE'])

    label = 'spatialDE' if ff == sde_files[-1] else None
    plt.plot(recall, precision, color=sde_color, label=label)

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.xlim(0, 1)
plt.ylim(0, 1)
axs[0].set_axisbelow(True)
plt.grid(color='#CCCCCC', dashes=[3, 3], lw=1)
plt.legend()

plt.sca(axs[1])
auprs = pd.DataFrame(auprs, columns=['AUPR', 'Method'])

sns.stripplot(x='Method', y='AUPR', data=auprs, jitter=True, palette=[hs_color, sde_color])
plt.ylim(0, .6)
plt.xlabel('')
axs[1].set_axisbelow(True)
plt.grid(color='#CCCCCC', axis='y', dashes=[3, 3], lw=1)
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.2)
# plt.show()
plt.savefig('AUPRs.svg')

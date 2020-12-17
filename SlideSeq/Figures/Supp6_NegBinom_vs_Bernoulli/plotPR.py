from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc
from functools import reduce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


hs_files = [
    "../../Puck_180819_9/hotspot/hotspot.txt",
    "../../Puck_180819_10/hotspot/hotspot.txt",
    "../../Puck_180819_11/hotspot/hotspot.txt",
    "../../Puck_180819_12/hotspot/hotspot.txt",
]

nb_files = [
    "../../Puck_180819_9/neg_binom/hotspot/hotspot.txt",
    "../../Puck_180819_10/neg_binom/hotspot/hotspot.txt",
    "../../Puck_180819_11/neg_binom/hotspot/hotspot.txt",
    "../../Puck_180819_12/neg_binom/hotspot/hotspot.txt",
]


# %% precision recall curve

hs_color = sns.color_palette("deep")[0]
sd_color = sns.color_palette("deep")[2]

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
    auprs.append([aupr, 'Hotspot-Bernoulli'])

    label = 'Hotspot-Bernoulli' if ff == hs_files[-1] else None
    plt.plot(recall, precision, color=hs_color, label=label)

for ff in nb_files:
    res = pd.read_table(ff, index_col=0)

    y_score = res.Z.values
    y_true = np.array(
        [x.upper() in all_markers_norm for x in res.index]
    ).astype('int')

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    aupr = auc(recall, precision)
    auprs.append([aupr, 'Hotspot-NegBinom'])

    label = 'Hotspot-NegBinom' if ff == nb_files[-1] else None
    plt.plot(recall, precision, color=sd_color, label=label)

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.xlim(0, 1)
plt.ylim(0, 1)
axs[0].set_axisbelow(True)
plt.grid(color='#CCCCCC', dashes=[3, 3], lw=1)
plt.legend()

plt.sca(axs[1])
auprs = pd.DataFrame(auprs, columns=['AUPR', 'Method'])

methods = ['Hotspot-Bernoulli', 'Hotspot-NegBinom']
sns.stripplot(
    x='Method', y='AUPR', data=auprs, jitter=True, palette=[hs_color, sd_color],
    order=methods
)
plt.ylim(0, .6)
plt.xlabel('')
plt.xticks([0, 1], [x.replace('-', '-\n') for x in methods])
axs[1].set_axisbelow(True)
plt.grid(color='#CCCCCC', axis='y', dashes=[3, 3], lw=1)
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.2)
plt.savefig('AUPRs.svg')

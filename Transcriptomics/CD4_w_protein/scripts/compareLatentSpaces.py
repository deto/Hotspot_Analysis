"""
Compare latent spaces using the protein autocorrelation
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%

l1 = pd.read_table("../hotspot/hotspot_threshold_proteins.txt", index_col=0)
l2 = pd.read_table("../hotspot/hotspot_hvg_proteins.txt", index_col=0)
l3 = pd.read_table("../hotspot/hotspot_hs_proteins.txt", index_col=0)
l4 = pd.read_table("../hotspot/hotspot_hvg2_proteins.txt", index_col=0)

l2 = l2.loc[l1.index]
l3 = l3.loc[l1.index]
l4 = l4.loc[l1.index]


plt.figure()
plt.hist(l1.C, 10)
plt.show()

plt.figure()
plt.plot(l1.Z, l1.C, 'o', ms=4)
plt.show()


delta = l2.C - l1.C
mean = (l2.C + l1.C)/2

plt.figure()
plt.plot(mean, delta, 'o', ms=4)
plt.show()

plt.figure()
plt.plot(l1.C, l2.C, 'o', ms=4)
plt.autoscale(False)
plt.plot([-100, 100], [-100, 100], '--', lw=1)
plt.show()

from bio_utils.plots import hover_plot

hover_plot(
    l1.C, l2.C, l1.index,
    'o', ms=4
)
plt.autoscale(False)
plt.plot([-100, 100], [-100, 100], '--', lw=1)
plt.show()

hover_plot(
    l2.C, l4.C, l2.index,
    'o', ms=4
)
plt.autoscale(False)
plt.plot([-100, 100], [-100, 100], '--', lw=1)
plt.show()


data = pd.concat(
    (
        l1.C.rename('Threshold'),
        l2.C.rename('Hvg'),
        l4.C.rename('Hvg2'),
        l3.C.rename('Hvg+Hs'),
    ), axis=1)
data.index.name = 'Protein'
data = data.loc[data.max(axis=1) > .05]
data = data.reset_index()
data['Hvg-Diff'] = (data['Hvg'] - data['Threshold']) / data['Threshold'] * 100
data['Hvg2-Diff'] = (data['Hvg2'] - data['Threshold']) / data['Threshold'] * 100
data['Hs-Diff'] = (data['Hvg+Hs'] - data['Threshold']) / data['Threshold'] * 100


data = data.loc[
    :, ['Protein', 'Hvg-Diff', 'Hvg2-Diff', 'Hs-Diff']
].melt(id_vars=['Protein'])


plt.figure()
sns.boxplot(data=data, x='variable', y='value')
plt.show()


# Great! Let's look at this across more datasets then
# %% Load dataset
ds = {}

ds["5K_All"] = {
    "thresh": pd.read_table(
        "../../PBMC_5k/hotspot/hotspot_threshold_proteins.txt", index_col=0
    ),
    "hvg": pd.read_table("../../PBMC_5k/hotspot/hotspot_hvg_proteins.txt", index_col=0),
    "hs": pd.read_table("../../PBMC_5k/hotspot/hotspot_hs_proteins.txt", index_col=0),
}

ds["5K_Mono"] = {
    "thresh": pd.read_table(
        "../../Mono_w_protein/hotspot/hotspot_threshold_proteins.txt", index_col=0
    ),
    "hvg": pd.read_table("../../Mono_w_protein/hotspot/hotspot_hvg_proteins.txt", index_col=0),
    "hs": pd.read_table("../../Mono_w_protein/hotspot/hotspot_hs_proteins.txt", index_col=0),
}

ds["5K_CD4"] = {
    "thresh": pd.read_table(
        "../hotspot/hotspot_threshold_proteins.txt", index_col=0
    ),
    "hvg": pd.read_table("../hotspot/hotspot_hvg_proteins.txt", index_col=0),
    "hs": pd.read_table("../hotspot/hotspot_hs_proteins.txt", index_col=0),
}

ds["Cite_CD4"] = {
    "thresh": pd.read_table(
        "../../CiteSeq_CD4/hotspot/hotspot_threshold_proteins.txt", index_col=0
    ),
    "hvg": pd.read_table("../../CiteSeq_CD4/hotspot/hotspot_hvg_proteins.txt", index_col=0),
    "hs": pd.read_table("../../CiteSeq_CD4/hotspot/hotspot_hs_proteins.txt", index_col=0),
}

ds["Cite_All"] = {
    "thresh": pd.read_table(
        "../../CiteSeq_All/hotspot/hotspot_threshold_proteins.txt", index_col=0
    ),
    "hvg": pd.read_table("../../CiteSeq_All/hotspot/hotspot_hvg_proteins.txt", index_col=0),
    "hs": pd.read_table("../../CiteSeq_All/hotspot/hotspot_hs_proteins.txt", index_col=0),
}

ds["10k_All"] = {
    "thresh": pd.read_table(
        "../../PBMC_10k/hotspot/hotspot_threshold_proteins.txt", index_col=0
    ),
    "hvg": pd.read_table("../../PBMC_10k/hotspot/hotspot_hvg_proteins.txt", index_col=0),
    "hs": pd.read_table("../../PBMC_10k/hotspot/hotspot_hs_proteins.txt", index_col=0),
}

ds["MALT"] = {
    "thresh": pd.read_table(
        "../../MALT/hotspot/hotspot_threshold_proteins.txt", index_col=0
    ),
    "hvg": pd.read_table("../../MALT/hotspot/hotspot_hvg_proteins.txt", index_col=0),
    "hs": pd.read_table("../../MALT/hotspot/hotspot_hs_proteins.txt", index_col=0),
}

# %% Plot results

def plot_dataset_xy(ds, y, x='thresh', ax=None):
    if ax is not None:
        plt.sca(ax)
    plt.plot(ds[x].C, ds[y].C, 'o', ms=4)
    plt.autoscale(False)
    plt.plot([-100, 100], [-100, 100], '--', lw=1)
    plt.ylim(-.1, 1.1)
    plt.xlim(-.1, 1.1)
    plt.xlabel(x)
    plt.ylabel(y)


def plot_dataset_box(ds, ax=None, cutoff=None):
    plt.sca(ax)
    ax = plt.gca()
    plt.grid(axis='y')
    data = pd.concat(
        (
            ds["thresh"].C.rename("Threshold"),
            ds["hvg"].C.rename("Hvg"),
            ds["hs"].C.rename("Hvg+Hs"),
        ),
        axis=1,
    )
    data.index.name = "Protein"

    if cutoff is not None:
        data = data.loc[data.max(axis=1) > cutoff]

    data = data.reset_index()
    data["Hvg-Diff"] = (
        (data["Hvg"] - data["Threshold"]) / data["Threshold"] * 100
    )
    data["Hs-Diff"] = (
        (data["Hvg+Hs"] - data["Threshold"]) / data["Threshold"] * 100
    )

    data = data.loc[:, ["Protein", "Hvg-Diff", "Hs-Diff"]].melt(
        id_vars=["Protein"]
    )

    sns.boxplot(data=data, x="variable", y="value", ax=ax)
    plt.ylim(-70, 70)
    plt.yticks(np.arange(-60, 61, 20))
    ax.set_axisbelow(True)


def dataset_metric(ds, cutoff=None):
    data = pd.concat(
        (
            ds["thresh"].C.rename("Threshold"),
            ds["hvg"].C.rename("Hvg"),
            ds["hs"].C.rename("Hvg+Hs"),
        ),
        axis=1,
    )
    data.index.name = "Protein"

    if cutoff is not None:
        data = data.loc[data.max(axis=1) > cutoff]

    data = data.reset_index()
    data["Hvg-Diff"] = (
        (data["Hvg"] - data["Threshold"]) / data["Threshold"] * 100
    )
    data["Hs-Diff"] = (
        (data["Hvg+Hs"] - data["Threshold"]) / data["Threshold"] * 100
    )

    return (data["Hvg-Diff"].mean(), data["Hs-Diff"].mean())



fig, axs = plt.subplots(len(ds), 4, figsize=(15, 9))

for ax_row, key in zip(axs, ds.keys()):
    plot_dataset_xy(ds[key], y='hvg', ax=ax_row[0])
    plot_dataset_xy(ds[key], y='hs', ax=ax_row[1])
    plot_dataset_box(ds[key], ax=ax_row[2])
    plot_dataset_box(ds[key], ax=ax_row[3], cutoff=0.05)
    plt.title(key)

    hvg, hs = dataset_metric(ds[key], cutoff=0.05)
    print("{}: Hvg: {:.2f}, HS: {:.2f}".format(key, hvg, hs))

plt.show()

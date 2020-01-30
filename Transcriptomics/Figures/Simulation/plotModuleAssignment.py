import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

plt.rcParams['svg.fonttype'] = 'none'

data_dir = "../../../data/Simulated4"
results_dir = "../../Simulation4"

def compute_stats(all_genes):

    f_confusion = (
        all_genes.groupby(["TrueModule", "PredictedModuleMapped"])
        .size()
        .rename("Size")
        .reset_index()
    )

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for row in f_confusion.itertuples():
        if row.TrueModule == -1:
            if row.PredictedModuleMapped == -1:
                tn += row.Size
            else:
                fp += row.Size
        else:
            if row.PredictedModuleMapped == row.TrueModule:
                tp += row.Size
            elif row.PredictedModuleMapped == -1:
                fn += row.Size
            else:
                fp += row.Size

    if tp+fp == 0:
        precision = float('nan')
    else:
        precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    accuracy = (tp + tn) / (tp + tn + fp + fn)

    return precision, recall, accuracy


# True Data Load data
def compute_results(rep, which):
    print(rep, which)

    rep_dir = "rep{}".format(rep)

    gene_indices = pd.read_table(
        os.path.join(data_dir, rep_dir, "gene_indices.txt"), index_col=0
    )
    gene_indices.head()

    gene_effects = pd.read_table(
        os.path.join(data_dir, rep_dir, "gene_effects.txt"), index_col=0
    )

    gene_means = pd.read_table(
        os.path.join(data_dir, rep_dir, "observed_counts.txt.gz"), index_col=0
    ).mean(axis=1)

    def gene_index_to_gene_id(ix):
        return "Gene{}".format(ix + 1)

    module_indices = {
        i + 1: gene_indices.iloc[:, i].map(gene_index_to_gene_id).values
        for i in range(gene_indices.shape[1])
    }

    module_indices = {i: pd.Index(v) for i, v in module_indices.items()}

    # Load predicted data
    if which == 'Hotspot':
        results = pd.read_table(
            os.path.join(results_dir, rep_dir, "hotspot/modules.txt"),
            index_col=0).Cluster
    elif which == 'Regular':
        results = pd.read_table(
            os.path.join(results_dir, rep_dir, "hotspot/modules_regular.txt"),
            index_col=0).Cluster
    else:
        raise Exception("Bad 'which' value!")

    # Combine
    all_genes = pd.DataFrame(index=gene_effects.index)
    all_genes["TrueModule"] = -1
    for m in module_indices:
        all_genes.loc[module_indices[m], "TrueModule"] = m

    all_genes["PredictedModule"] = -1
    for m in results.unique():
        all_genes.loc[results.index[results == m], "PredictedModule"] = m

    # Optional - do we subset to only the genes already selected for 'pairs'?
    all_genes = all_genes.loc[results.index]

    # Now, need to map indices between the two
    # For each true module, gets top two 'predicted modules'
    # Raise an error if we are re-assigning something that has been assigned
    pred_to_true = {}
    pred_to_true[-1] = -1
    true_to_pred = {}

    confusion = (
        all_genes.groupby(["TrueModule", "PredictedModule"]).size().reset_index()
    )

    confusion = confusion.loc[confusion.PredictedModule != -1]
    confusion = confusion.sort_values(0, ascending=False)

    for _, row in confusion.iterrows():
        if row.TrueModule not in true_to_pred:
            true_to_pred[row.TrueModule] = []

        if len(true_to_pred[row.TrueModule]) == 2:
            continue

        if row.PredictedModule in pred_to_true:
            continue

        pred_to_true[row.PredictedModule] = row.TrueModule
        true_to_pred[row.TrueModule].append(row.PredictedModule)

    all_genes["PredictedModuleMapped"] = all_genes.PredictedModule.map(pred_to_true)

    precision, recall, accuracy = compute_stats(all_genes)

    # What about for gene_effect bins?
    ge_sub = gene_effects.loc[results.index]
    ge_size = ge_sub.iloc[:, 1:].abs().max(axis=1)
    ge_size = ge_size.loc[ge_size > 0]
    ge_groups = pd.qcut(ge_size, 5)

    group_results = []
    for i, cc in enumerate(ge_groups.cat.categories):
        genes = ge_groups.index[ge_groups == cc]
        p, r, a = compute_stats(all_genes.loc[genes])
        group_results.append(
            [i, p, r, a]
        )
    group_results = pd.DataFrame(
        group_results,
        columns=['GEQ', 'Precision', 'Recall', 'Accuracy']
    )

    # What about for each module?
    ge_sub = gene_effects.loc[results.index]

    module_results = []
    for i, mod in enumerate(all_genes.TrueModule.unique()):
        if mod == -1:
            continue
        p, r, a = compute_stats(all_genes.loc[all_genes.TrueModule == mod])
        module_results.append(
            [mod, p, r, a]
        )
    module_results = pd.DataFrame(
        module_results,
        columns=['TrueModule', 'Precision', 'Recall', 'Accuracy']
    )

    # What about for gene_mean bins?
    ge_sub = gene_means.loc[results.index]
    ge_groups = pd.qcut(ge_sub, 5)

    mean_results = []
    for i, cc in enumerate(ge_groups.cat.categories):
        genes = ge_groups.index[ge_groups == cc]
        p, r, a = compute_stats(all_genes.loc[genes])
        mean_results.append(
            [i, p, r, a]
        )
    mean_results = pd.DataFrame(
        mean_results,
        columns=['Mean', 'Precision', 'Recall', 'Accuracy']
    )

    # What about combination of mean and of gene effect size?
    ge_sub = gene_effects.loc[results.index]
    ge_size = ge_sub.iloc[:, 1:].abs().max(axis=1)
    ge_size = ge_size.loc[ge_size > 0]
    ge_groups_ge = pd.qcut(ge_size, 5)
    ge_sub = gene_means.loc[results.index]
    ge_groups_mu = pd.qcut(ge_sub, 5)

    ge_mean_results = []
    for i, cc_ge in enumerate(ge_groups_ge.cat.categories):
        for j, cc_mu in enumerate(ge_groups_mu.cat.categories):
            genes = (
                (ge_groups_ge.index[ge_groups_ge == cc_ge]) &
                (ge_groups_mu.index[ge_groups_mu == cc_mu])
            )

            p, r, a = compute_stats(all_genes.loc[genes])

            ge_mean_results.append(
                [rep, i, j, p, r, a, genes.size]
            )

    ge_mean_results = pd.DataFrame(
        ge_mean_results,
        columns=['Rep', 'GEQ', 'MeanQuantile', 'Precision', 'Recall', 'Accuracy', 'NGenes']
    )

    return (precision, recall, accuracy), group_results, module_results, mean_results, ge_mean_results

# %%

reps = list(range(1, 11))
from tqdm import tqdm

all_res_hs = [compute_results(x, 'Hotspot') for x in tqdm(reps)]
group_results_hs = [x[1] for x in all_res_hs]
module_results_hs = [x[2] for x in all_res_hs]
mean_results_hs = [x[3] for x in all_res_hs]
ge_mean_results_hs = [x[4] for x in all_res_hs]
all_res_hs = [x[0] for x in all_res_hs]
all_res_hs = pd.DataFrame(
    all_res_hs, columns=['Precision', 'Recall', 'Accuracy']
)
all_res_hs['Method'] = 'Hotspot'

group_results_hs = pd.concat(
    group_results_hs, axis=0
)
group_results_hs['Method'] = 'Hotspot'

module_results_hs = pd.concat(
    module_results_hs, axis=0
)
module_results_hs['Method'] = 'Hotspot'

mean_results_hs = pd.concat(
    mean_results_hs, axis=0
)
mean_results_hs['Method'] = 'Hotspot'

ge_mean_results_hs = pd.concat(
    ge_mean_results_hs, axis=0
)
ge_mean_results_hs['Method'] = 'Hotspot'

all_res_reg = [compute_results(x, 'Regular') for x in tqdm(reps)]
group_results_reg = [x[1] for x in all_res_reg]
module_results_reg = [x[2] for x in all_res_reg]
mean_results_reg = [x[3] for x in all_res_reg]
ge_mean_results_reg = [x[4] for x in all_res_reg]
all_res_reg = [x[0] for x in all_res_reg]
all_res_reg = pd.DataFrame(
    all_res_reg, columns=['Precision', 'Recall', 'Accuracy']
)
all_res_reg['Method'] = 'Regular'

group_results_reg = pd.concat(
    group_results_reg, axis=0
)
group_results_reg['Method'] = 'Regular'

module_results_reg = pd.concat(
    module_results_reg, axis=0
)
module_results_reg['Method'] = 'Regular'

mean_results_reg = pd.concat(
    mean_results_reg, axis=0
)
mean_results_reg['Method'] = 'Regular'

ge_mean_results_reg = pd.concat(
    ge_mean_results_reg, axis=0
)
ge_mean_results_reg['Method'] = 'Regular'

all_res = pd.concat(
    (all_res_hs, all_res_reg),
)

group_res = pd.concat(
    (group_results_hs, group_results_reg),
)

module_res = pd.concat(
    (module_results_hs, module_results_reg),
)

mean_res = pd.concat(
    (mean_results_hs, mean_results_reg),
)

ge_mean_res = pd.concat(
    (ge_mean_results_hs, ge_mean_results_reg),
)

# %%


# plt.figure()
# sns.boxplot(data=all_res, x='Method', y='Precision')
# plt.show()
# 
# plt.figure()
# sns.boxplot(data=all_res, x='Method', y='Recall')
# plt.show()
# 
# plt.figure()
# sns.boxplot(data=all_res, x='Method', y='Accuracy')
# plt.show()
# 
# 
# plt.figure()
# sns.boxplot(data=group_res, x='GEQ', y='Accuracy', hue='Method', whis='range')
# plt.show()
# 
# plt.figure()
# sns.boxplot(data=group_res, x='GEQ', y='Recall', hue='Method')
# plt.show()


# %%

fig, axs = plt.subplots(1, 2, gridspec_kw=dict(width_ratios=[0.5, 1]), sharey=True)

m_order = ['Regular', 'Hotspot']

plt.sca(axs[0])
axs[0].set_axisbelow(True)
sns.stripplot(data=all_res, x='Method', y='Accuracy', order=m_order)
plt.title('Overall\n Accuracy')
plt.grid(color='#999999', ls='--', lw=.5, axis='y')
plt.sca(axs[1])
axs[1].set_axisbelow(True)
sns.barplot(
    data=group_res, x='GEQ', y='Accuracy',
    hue='Method', hue_order=m_order
)
plt.xlabel('Quantile\n(Low Effect Size to High)')
plt.xticks([0, 1, 2, 3, 4], ['1', '2', '3', '4', '5'])
plt.title('Accuracy per\n Gene Effect Quantile')
plt.grid(axis='y', color='#999999', ls='--', lw=.5)
plt.subplots_adjust(bottom=0.2)
# plt.show()
plt.savefig('ModuleAccuracy.svg')

# %% By Module

# plt.figure()
# sns.barplot(
#     data=module_res, x='TrueModule', y='Accuracy',
#     hue='Method', hue_order=m_order
# )
# plt.xlabel('Gene Module')
# #plt.xticks([0, 1, 2, 3, 4], ['1', '2', '3', '4', '5'])
# plt.title('Accuracy per\n Gene Effect Quantile')
# plt.grid(axis='y', color='#999999', ls='--', lw=.5)
# plt.grid(axis='x', color='#999999', ls='--', lw=.5, which='minor')
# plt.subplots_adjust(bottom=0.2)
# plt.show()

# %% By Mean

# mean_res.columns = ['MeanQuantile', 'Precision', 'Recall', 'Accuracy', 'Method']
# 
# plt.figure()
# sns.barplot(
#     data=mean_res, x='MeanQuantile', y='Accuracy',
#     hue='Method', hue_order=m_order
# )
# plt.xlabel('Mean Quantile')
# #plt.xticks([0, 1, 2, 3, 4], ['1', '2', '3', '4', '5'])
# plt.title('Accuracy per\n Gene Mean Quantile')
# plt.grid(axis='y', color='#999999', ls='--', lw=.5)
# plt.grid(axis='x', color='#999999', ls='--', lw=.5, which='minor')
# plt.subplots_adjust(bottom=0.2)
# plt.show()


# %% By combinations of mean and effect size:

if ge_mean_res['GEQ'].min() == 0:
    ge_mean_res['GEQ'] = ge_mean_res['GEQ'] + 1

if ge_mean_res['MeanQuantile'].min() == 0:
    ge_mean_res['MeanQuantile'] = ge_mean_res['MeanQuantile'] + 1

plot_data_hs = ge_mean_res.loc[ge_mean_res.Method == 'Hotspot']
plot_data_hs = plot_data_hs.groupby(['GEQ', 'MeanQuantile']).mean() \
    .reset_index() \
    .pivot(index='GEQ', columns='MeanQuantile', values='Accuracy')

plot_data_reg = ge_mean_res.loc[ge_mean_res.Method == 'Regular']
plot_data_reg = plot_data_reg.groupby(['GEQ', 'MeanQuantile']).mean() \
    .reset_index() \
    .pivot(index='GEQ', columns='MeanQuantile', values='Accuracy')


delta = (
    ge_mean_res.loc[ge_mean_res.Method == 'Hotspot']
    .set_index(['Rep', 'GEQ', 'MeanQuantile'])['Accuracy']
) - (
    ge_mean_res.loc[ge_mean_res.Method == 'Regular']
    .set_index(['Rep', 'GEQ', 'MeanQuantile'])['Accuracy']
)

delta = delta.reset_index().groupby(['GEQ', 'MeanQuantile'])['Accuracy'].mean()
delta = delta.reset_index().pivot(
    index='GEQ', columns='MeanQuantile', values='Accuracy'
)

fig, axs = plt.subplots(2, 2, figsize=(7, 7),
                        gridspec_kw=dict(
                            hspace=0.55, wspace=0.35
                        )
                        )
plt.sca(axs[0, 0])
sns.heatmap(plot_data_hs*100, vmin=0, vmax=100,
            cbar_kws=dict(
                ticks=[0, 50, 100],
            ))
plt.title('Hotspot')
plt.xlabel('Gene Expression Level\n(Quantile)')
plt.ylabel('Gene-Effect\n(Quantile)')

plt.sca(axs[0, 1])
sns.heatmap(plot_data_reg*100, vmin=0, vmax=100,
            cbar_kws=dict(
                ticks=[0, 50, 100],
            ))
plt.title('Pearson')
plt.xlabel('Gene Expression Level\n(Quantile)')
plt.ylabel('Gene-Effect\n(Quantile)')

plt.sca(axs[1, 0])
sns.heatmap(delta*100, cmap="RdBu_r",
            vmin=-30, vmax=30,
            cbar_kws=dict(
                ticks=[-30, 0, 30],
            ))
plt.xlabel('Gene Expression Level\n(Quantile)')
plt.ylabel('Gene-Effect\n(Quantile)')

axs[1, 1].remove()
plt.savefig('ModuleAccuracy2.svg')

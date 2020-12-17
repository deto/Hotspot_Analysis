import os
from sklearn.metrics import roc_curve, precision_recall_curve, average_precision_score, roc_auc_score
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
import seaborn as sns
import matplotlib.ticker as ticker

plt.rcParams['svg.fonttype'] = 'none'

results_dir = "../../Simulation6"
data_dir = "../../../data/Simulated6"

def get_module_indices(directory):
    """
    Returns the gene-module assignment mapping
    """

    try:

        gene_indices = pd.read_table(
            os.path.join(directory, "gene_indices.txt"), index_col=0
        )
        gene_indices.head()

        def gene_index_to_gene_id(ix):
            return "Gene{}".format(ix + 1)

        module_indices = {
            i + 1: gene_indices.iloc[:, i].map(gene_index_to_gene_id).values
            for i in range(gene_indices.shape[1])
        }

        module_indices = {i: pd.Index(v) for i, v in module_indices.items()}

    except OSError:

        gene_effects = pd.read_table(
            os.path.join(directory, "gene_effects.txt"), index_col=0
        )

        gene_effects = gene_effects.iloc[:, 1:]

        module_indices = {}
        for i in range(gene_effects.shape[1]):
            vals = gene_effects.iloc[:, i]
            module_indices[i+1] = vals[vals != 0].index

    return module_indices


def get_results_hs(rep, k):

    rep_dir = "rep{}/".format(rep)

    hs_results_file = os.path.join(
        results_dir, rep_dir, "knn_sensitivity/k_{k}/hotspot/hotspot.txt".format(k=k)
    )

    hs_results = pd.read_table(hs_results_file, index_col=0)

    module_indices = get_module_indices(os.path.join(data_dir, rep_dir))

    all_module_indices = set()
    for x in module_indices.values():
        all_module_indices |= set(x)
    all_module_indices = pd.Index(all_module_indices)

    # Plot an ROC curve
    hs_results["True Positive"] = [
        1 if x in all_module_indices else 0 for x in hs_results.index
    ]

    fpr, tpr, thresholds = roc_curve(
        hs_results["True Positive"], hs_results["Z"]
    )

    precision, recall, thresholds = precision_recall_curve(
        hs_results["True Positive"], hs_results["Z"]
    )

    auprc = average_precision_score(
        hs_results["True Positive"], hs_results["Z"]
    )

    auroc = roc_auc_score(
        hs_results["True Positive"], hs_results["Z"]
    )

    # Also, use the FDR to compute point-wise estimate
    pos = hs_results["True Positive"].sum()
    neg = hs_results.shape[0] - pos
    hs_sub = hs_results.loc[hs_results.FDR < .05]
    tp = hs_sub["True Positive"].sum()
    fp = hs_sub.shape[0] - tp
    fn = pos - tp
    tn = hs_results.shape[0] - tp - fp - fn

    point_accuracy = (tp + tn) / (tp + tn + fp + fn)
    point_precision = tp / (tp + fp)
    point_recall = tp / (tp + fn)
    point_fpr = fp / neg

    out = {
        'fpr': fpr,
        'tpr': tpr,
        'precision': precision,
        'recall': recall,
        'point_est': {
            'accuracy': point_accuracy,
            'precision': point_precision,
            'recall': point_recall,
            'tpr': point_recall,
            'fpr': point_fpr,
            'k': k,
        },
        'summary_stats': {
            'auroc': auroc,
            'auprc': auprc,
            'k': k,
        },
        'k': k
    }

    return out


# %%

rep = 1
k_values = [5, 10, 30, 50, 100, 300, 500, 1000]
res_hs = {k: get_results_hs(rep, k) for k in k_values}

point_estimates = [pd.Series(x['point_est']) for x in res_hs.values()]
point_estimates = pd.concat(point_estimates, axis=1).T

summary_stats = [pd.Series(x['summary_stats']) for x in res_hs.values()]
summary_stats = pd.concat(summary_stats, axis=1).T

# %% Plot AUROC

plt.figure()
ax = plt.gca()

plt.plot(
    summary_stats.k, summary_stats.auroc, 'o-', label='3000 cells'
)

plt.xlabel('# of Neighbors')
plt.ylabel('Area Under ROC Curve')
plt.xscale('log')
plt.xticks(summary_stats.k)

ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_ticks([], minor=True)
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', ls=(0, (5, 5)), lw=.5)
plt.savefig('AU_ROC_k_sensitivity.svg', dpi=300)


# %% Plot AUPRC

plt.figure()
ax = plt.gca()

plt.plot(
    summary_stats.k, summary_stats.auprc, 'o-', label='3000 cells'
)

plt.xlabel('# of Neighbors')
plt.ylabel('Area Under\nPrecision-Recall Curve')
plt.xscale('log')
plt.xticks(summary_stats.k)

ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_ticks([], minor=True)
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', ls=(0, (5, 5)), lw=.5)
plt.savefig('AU_PR_k_sensitivity.svg', dpi=300)

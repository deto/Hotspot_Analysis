import os
from sklearn.metrics import roc_curve, precision_recall_curve
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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


def get_results_hs(rep):

    rep_dir = "rep{}/".format(rep)

    hs_results = pd.read_table(
        os.path.join(results_dir, rep_dir, "hotspot/hotspot_threshold.txt"),
        index_col=0
    )

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
        }
    }

    return out


def get_results_hvg(rep):
    rep_dir = "rep{}/".format(rep)

    hvg_results = pd.read_table(
        os.path.join(results_dir, rep_dir, "genes/hvg_info.txt"),
        index_col=0
    )

    module_indices = get_module_indices(os.path.join(data_dir, rep_dir))

    all_module_indices = set()
    for x in module_indices.values():
        all_module_indices |= set(x)
    all_module_indices = pd.Index(all_module_indices)

    # Plot an ROC curve
    hvg_results["True Positive"] = [
        1 if x in all_module_indices else 0 for x in hvg_results.index
    ]

    fpr, tpr, thresholds = roc_curve(
        hvg_results["True Positive"], hvg_results["gene.dispersion.scaled"]
    )

    precision, recall, thresholds = precision_recall_curve(
        hvg_results["True Positive"], hvg_results["gene.dispersion.scaled"]
    )

    out = {
        'fpr': fpr,
        'tpr': tpr,
        'precision': precision,
        'recall': recall,
    }

    return out


def get_results_danb(rep):
    rep_dir = "rep{}/".format(rep)

    danb_results = pd.read_table(
        os.path.join(results_dir, rep_dir, "genes/danb_info.txt"),
        index_col=0
    )

    module_indices = get_module_indices(os.path.join(data_dir, rep_dir))

    all_module_indices = set()
    for x in module_indices.values():
        all_module_indices |= set(x)
    all_module_indices = pd.Index(all_module_indices)

    # Plot an ROC curve
    danb_results["True Positive"] = [
        1 if x in all_module_indices else 0 for x in danb_results.index
    ]

    fpr, tpr, thresholds = roc_curve(
        danb_results["True Positive"], danb_results["p.value"]*-1
    )

    precision, recall, thresholds = precision_recall_curve(
        danb_results["True Positive"], danb_results["p.value"]*-1
    )

    out = {
        'fpr': fpr,
        'tpr': tpr,
        'precision': precision,
        'recall': recall,
    }

    return out


def get_results_pca(rep):
    rep_dir = "rep{}/".format(rep)

    pca_results = pd.read_table(
        os.path.join(results_dir, rep_dir, "genes/pca_info.txt"),
        index_col=0
    )

    module_indices = get_module_indices(os.path.join(data_dir, rep_dir))

    all_module_indices = set()
    for x in module_indices.values():
        all_module_indices |= set(x)
    all_module_indices = pd.Index(all_module_indices)

    # Plot an ROC curve
    pca_results["True Positive"] = [
        1 if x in all_module_indices else 0 for x in pca_results.index
    ]

    fpr, tpr, thresholds = roc_curve(
        pca_results["True Positive"], pca_results["Score"]
    )

    precision, recall, thresholds = precision_recall_curve(
        pca_results["True Positive"], pca_results["Score"]
    )

    out = {
        'fpr': fpr,
        'tpr': tpr,
        'precision': precision,
        'recall': recall,
    }

    return out


def gather_results(reps, fxn):

    all_results_pre = [
        fxn(i) for i in reps
    ]

    # Gather tpr
    all_tpr = []
    new_fpr = np.linspace(0, 1, 500)
    for rr in all_results_pre:
        tpr = rr["tpr"]
        fpr = rr["fpr"]

        # need to interpolate tpr as function of fpr so we can combine
        new_tpr = np.interp(x=new_fpr, xp=fpr, fp=tpr)

        # You can have multiple tpr values at fpr=0
        # Just set tpr[fpr=0] = 0 for plotting sake
        new_tpr[0] = 0

        all_tpr.append(new_tpr)
    all_tpr = np.vstack(all_tpr)
    tpr_means = all_tpr.mean(axis=0)
    tpr_sds = all_tpr.std(axis=0)
    tpr_maxs = all_tpr.max(axis=0)
    tpr_mins = all_tpr.min(axis=0)

    # Gather precision
    all_pre = []
    new_recall = np.linspace(0, 1, 500)
    for rr in all_results_pre:
        precision = rr["precision"][::-1]
        recall = rr["recall"][::-1]

        # need to interpolate tpr as function of fpr so we can combine
        new_pre = np.interp(x=new_recall, xp=recall, fp=precision)

        all_pre.append(new_pre)

    all_pre = np.vstack(all_pre)
    pre_means = all_pre.mean(axis=0)
    pre_sds = all_pre.std(axis=0)
    pre_maxs = all_pre.max(axis=0)
    pre_mins = all_pre.min(axis=0)

    # Gather point est
    point_precision = []
    point_recall = []
    for rr in all_results_pre:
        if "point_est" in rr:
            precision = rr["point_est"]["precision"]
            recall = rr["point_est"]["recall"]
        else:
            precision = None
            recall = None

        point_precision.append(precision)
        point_recall.append(recall)

    out = {
        'all_res': all_results_pre,
        'fpr': new_fpr,
        'all_tpr': all_tpr,
        'tpr_means': tpr_means,
        'tpr_sds': tpr_sds,
        'tpr_maxs': tpr_maxs,
        'tpr_mins': tpr_mins,
        'recall': new_recall,
        'pre_means': pre_means,
        'pre_sds': pre_sds,
        'pre_maxs': pre_maxs,
        'pre_mins': pre_mins,
        'point_precision': point_precision,
        'point_recall': point_recall,
    }

    return out


# %%

reps = [x+1 for x in range(10)]
res_hs = gather_results(reps, get_results_hs)
res_hvg = gather_results(reps, get_results_hvg)
res_danb = gather_results(reps, get_results_danb)
res_pca = gather_results(reps, get_results_pca)

# %% AUC plots

plt.figure()

fpr = res_hs['fpr']
mean = res_hs['tpr_means']
maxs = res_hs['tpr_maxs']
mins = res_hs['tpr_mins']


plt.plot(fpr, mean, "-", label='HS')
plt.fill_between(fpr, mins, maxs, color="gray", alpha=0.2)

fpr = res_hvg['fpr']
mean = res_hvg['tpr_means']
maxs = res_hvg['tpr_maxs']
mins = res_hvg['tpr_mins']


plt.plot(fpr, mean, "-", label='Hvg')
plt.fill_between(fpr, mins, maxs, color="gray", alpha=0.2)

fpr = res_danb['fpr']
mean = res_danb['tpr_means']
maxs = res_danb['tpr_maxs']
mins = res_danb['tpr_mins']


plt.plot(fpr, mean, "-", label='DANB')
plt.fill_between(fpr, mins, maxs, color="gray", alpha=0.2)

fpr = res_pca['fpr']
mean = res_pca['tpr_means']
maxs = res_pca['tpr_maxs']
mins = res_pca['tpr_mins']


plt.plot(fpr, mean, "-", label='PCA')
plt.fill_between(fpr, mins, maxs, color="gray", alpha=0.2)

plt.legend()
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.title('ROC Curve')
# plt.show()
plt.savefig('AUC.svg')

# %% PR plots

plt.figure()

colors = sns.color_palette('deep')

recall = res_hs['recall']
mean = res_hs['pre_means']
mins = res_hs['pre_mins']
maxs = res_hs['pre_maxs']

plt.plot(recall, mean, "-", label='Hotspot', color=colors[0])
plt.fill_between(recall, mins, maxs, color="gray", alpha=0.2, linewidth=0)

recall = res_hvg['recall']
mean = res_hvg['pre_means']
mins = res_hvg['pre_mins']
maxs = res_hvg['pre_maxs']

plt.plot(recall, mean, "-", label='HVG', color=colors[1])
plt.fill_between(recall, mins, maxs, color="gray", alpha=0.2, linewidth=0)

recall = res_danb['recall']
mean = res_danb['pre_means']
mins = res_danb['pre_mins']
maxs = res_danb['pre_maxs']

plt.plot(recall, mean, "-", label='NBDisp', color=colors[2])
plt.fill_between(recall, mins, maxs, color="gray", alpha=0.2, linewidth=0)

recall = res_pca['recall']
mean = res_pca['pre_means']
mins = res_pca['pre_mins']
maxs = res_pca['pre_maxs']

plt.plot(recall, mean, "-", label='PCA', color=colors[3])
plt.fill_between(recall, mins, maxs, color="gray", alpha=0.2, linewidth=0)

plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curves')

plt.grid(color='#dddddd', linestyle=(0, (5, 5)))
plt.gca().set_axisbelow(True)
# plt.show()
plt.savefig('PR.svg')

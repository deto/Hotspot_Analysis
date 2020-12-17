import os
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

base_dir = '../../Mono_w_protein'

datasets = [
    {
        'Name': 'Hotspot',
        'Train': os.path.join(base_dir, 'train/hotspot_hs/modules.txt'),
        'Test': os.path.join(base_dir, 'test/hotspot_hs/modules.txt'),
    },
    {
        'Name': 'WGCNA',
        'Train': os.path.join(base_dir, 'train/wgcna_hs/modules.txt'),
        'Test': os.path.join(base_dir, 'test/wgcna_hs/modules.txt'),
    },
    {
        'Name': 'ICA5',
        'Train': os.path.join(base_dir, 'train/ica5/modules.txt'),
        'Test': os.path.join(base_dir, 'test/ica5/modules.txt'),
    },
    # # Don't need so many ICA versions as this one isn't that different from ICA10
    # {
    #     'Name': 'ICA8',
    #     'Train': os.path.join(base_dir, 'train/ica8/modules.txt'),
    #     'Test': os.path.join(base_dir, 'test/ica8/modules.txt'),
    # },
    {
        'Name': 'ICA10',
        'Train': os.path.join(base_dir, 'train/ica10/modules.txt'),
        'Test': os.path.join(base_dir, 'test/ica10/modules.txt'),
    },
    {
        'Name': 'Grnboost',
        'Train': os.path.join(base_dir, 'train/grnboost/modules.txt'),
        'Test': os.path.join(base_dir, 'test/grnboost/modules.txt'),
    },
]

for data in datasets:
    train = pd.read_table(data['Train'], index_col=0)
    test = pd.read_table(data['Test'], index_col=0)

    train = train.Cluster.to_dict()
    test = test.Cluster.to_dict()
    data['TrainDict'] = train
    data['TestDict'] = test


def eval_module_consistency(data):
    consistency = (
        eval_module_consistency_inner(data['TrainDict'], data['TestDict']) +
        eval_module_consistency_inner(data['TestDict'], data['TrainDict'])
    ) / 2

    data['Consistency'] = consistency


def eval_module_consistency_inner(dict_a, dict_b):

    all_genes = set(dict_a.keys()) & set(dict_b.keys())

    # For each pairs of genes that are in the same module in 'A', how many are in the same module in 'B'?

    denom = 0
    num = 0
    for ga in all_genes:
        for gb in all_genes:

            if ga == gb: continue

            if dict_a[ga] == dict_a[gb] and dict_a[ga] != -1 and dict_a[gb] != -1:  # Same module in A
                denom += 1

                if dict_b[ga] == dict_b[gb] and dict_b[ga] != -1 and dict_b[gb] != -1:  # Same module in B
                    num += 1

    num = num/2
    denom = denom/2

    consistent_rate = num/denom

    return consistent_rate


def eval_num_modules(data):
    num_modules = (
        pd.Series(data['TrainDict']).unique().size - 1 +
        pd.Series(data['TestDict']).unique().size - 1
    ) / 2

    data['NumModules'] = num_modules


def eval_num_assigned(data):
    assigned = (
        (pd.Series(data['TrainDict']) != -1).sum()/2 +
        (pd.Series(data['TestDict']) != -1).sum()/2
    )

    data['NumAssigned'] = assigned


for data in tqdm(datasets):
    train = pd.read_table(data['Train'], index_col=0)
    test = pd.read_table(data['Test'], index_col=0)

    train = train.Cluster.to_dict()
    test = test.Cluster.to_dict()
    data['TrainDict'] = train
    data['TestDict'] = test

    eval_module_consistency(data)
    eval_num_modules(data)
    eval_num_assigned(data)


# %% Consolidate into a nice dataframe
columns = [
    'Name',
    'Consistency',
    'NumModules',
    'NumAssigned'
]

results = []
for data in datasets:
    results.append(
        [data[x] for x in columns]
    )

results = pd.DataFrame(results, columns=columns)


# %% Plot it
order = ['ICA5', 'ICA10', 'Grnboost', 'WGCNA', 'Hotspot']
colors = sns.color_palette('deep')[:len(order)]
plot_data = results.set_index('Name').loc[order]

fig, axs = plt.subplots(1, 3, figsize=(12, 4))

plt.sca(axs[0])

plt.bar(
    plot_data.index, plot_data.Consistency, alpha=0.9, color=colors
)
plt.xticks(rotation=45)
plt.ylabel('Proportion of Gene Pairs Which\nReplicate Across Data Split')
plt.title('Reproducibility')
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', lw=0.5, axis='y', ls=(0, (5, 5)))

plt.sca(axs[1])

plt.bar(
    plot_data.index, plot_data.NumModules, alpha=0.9, color=colors
)
plt.xticks(rotation=45)
plt.ylabel('Modules')
plt.title('# Modules')
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', lw=0.5, axis='y', ls=(0, (5, 5)))

plt.sca(axs[2])

plt.bar(
    plot_data.index, plot_data.NumAssigned, alpha=0.9, color=colors
)
plt.xticks(rotation=45)
plt.ylabel('Genes')
plt.title('# Genes Assigned')
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', lw=0.5, axis='y', ls=(0, (5, 5)))

plt.subplots_adjust(bottom=.25, wspace=0.4, left=0.1, right=0.9)
plt.savefig('Monocyte_Module_TrainTest.svg')

import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %% Need to gather all the timings

data = []
for results in glob.glob("../../Puck_180819_12/hotspot/timing/*.txt"):
    data.append(pd.read_table(results))

data = pd.concat(data, axis=0)

data_sd = []
for results in glob.glob("../../Puck_180819_12/spatialDE/timing/*.txt"):
    data_sd.append(pd.read_table(results))

data_sd = pd.concat(data_sd, axis=0)

# %%

plt.figure()
ax = plt.gca()
for genes, group in data.groupby("Genes"):
    group = group.sort_values('Cells')
    plt.plot(group.Cells, group.ElapsedSeconds, 'o-', label=genes)

#ax.set_xscale('log')
ax.set_xticks([5000, 10000, 20000])
ax.set_xticklabels(["5000", "10000", "20000"])
plt.legend()
plt.setp(ax.get_xminorticklabels(), visible=False)
#ax.get_xaxis().set_tick_params(which='minor', size=0)
#ax.get_xaxis().set_tick_params(which='minor', width=0)
plt.show()

# %%

plt.figure()
ax = plt.gca()
for genes, group in data_sd.groupby("Genes"):
    group = group.sort_values('Cells')
    plt.plot(group.Cells, group.SDESeconds, 'o-', label=genes)

#ax.set_xscale('log')
ax.set_xticks([5000, 10000, 20000])
ax.set_xticklabels(["5000", "10000", "20000"])
plt.legend()
plt.setp(ax.get_xminorticklabels(), visible=False)
#ax.get_xaxis().set_tick_params(which='minor', size=0)
#ax.get_xaxis().set_tick_params(which='minor', width=0)
plt.show()

# %% Plot them together

cmap = plt.get_cmap('tab10').colors
gene_cmap = {
    500: cmap[0],
    1000: cmap[1],
    2000: cmap[2],
    4000: cmap[3],
}

plt.figure()
ax = plt.gca()
lines = []
for genes, group in data.groupby("Genes"):
    group = group.sort_values('Cells')
    l = plt.plot(group.Cells, group.ElapsedSeconds, 'o-',
             label=genes, color=gene_cmap[genes])[0]
    lines.append(l)

for genes, group in data_sd.groupby("Genes"):
    group = group.sort_values('Cells')
    plt.plot(group.Cells, group.SDESeconds, 'o--',
             label=genes, color=gene_cmap[genes])

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xticks([5000, 10000, 20000])
ax.set_xticklabels(["5000", "10000", "20000"])
plt.setp(ax.get_xminorticklabels(), visible=False)
#ax.get_xaxis().set_tick_params(which='minor', size=0)
#ax.get_xaxis().set_tick_params(which='minor', width=0)
plt.grid(color='#cccccc')
plt.xlabel('# of Cells')
plt.ylabel('Runtime (s)')
plt.subplots_adjust(right=.78)

from matplotlib.legend import Legend
from matplotlib.lines import Line2D
lines = lines[::-1]
ll = Legend(ax, lines, [l.get_label() for l in lines], loc='upper left',
            bbox_to_anchor=(1, .8), title='# Genes')
ax.add_artist(ll)

lines = [
    Line2D([], [], linestyle='solid', label='Hotspot', color='#555555'),
    Line2D([], [], linestyle='dashed', label='SpatialDE', color='#555555'),
]
lines = lines[::-1]
ll = Legend(ax, lines, [l.get_label() for l in lines], loc='upper left',
            bbox_to_anchor=(1, 1), title='Method')
ax.add_artist(ll)

plt.show()
#plt.savefig('Timing.svg')

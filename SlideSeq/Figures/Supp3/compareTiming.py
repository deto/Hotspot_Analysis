import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json

plt.rcParams['svg.fonttype'] = 'none'

# %%

# %% Need to gather all the timings

data = []
for results in glob.glob("../../Puck_180819_12/hotspot/timing_pairs/*.txt"):
    data.append(pd.read_table(results))

data = pd.concat(data, axis=0)

data = data.loc[data.Genes == 500]

data_sd = []
for results in glob.glob("../../Puck_180819_12/spatialDE/timing_pairs/*.json"):
    xx = json.load(open(results))
    xx = pd.Series(xx).to_frame().T
    data_sd.append(xx)

data_sd = pd.concat(data_sd, axis=0)

# %% Plot them together

colors = plt.get_cmap('tab10').colors

data = data.sort_values('ElapsedSeconds')
data_sd = data_sd.sort_values('Time')

plt.figure()
ax = plt.gca()
plt.plot(data.Cells, data.ElapsedSeconds, "o-", label="Hotspot", color=colors[0])

plt.plot(data_sd.Cells, data_sd.Time, "o--", label="SpatialDE", color=colors[0])

ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xticks([1000, 2000, 5000, 10000])
ax.set_xticklabels(["1000", "2000", "5000", "10000"])
plt.setp(ax.get_xminorticklabels(), visible=False)
# ax.get_xaxis().set_tick_params(which='minor', size=0)
# ax.get_xaxis().set_tick_params(which='minor', width=0)
plt.grid(color="#cccccc")
plt.xlabel("# of Cells")
plt.ylabel("Runtime (s)")
plt.subplots_adjust(right=0.78)

from matplotlib.legend import Legend
from matplotlib.lines import Line2D
lines = [
    Line2D([], [], linestyle='solid', label='Hotspot', color=colors[0]),
    Line2D([], [], linestyle='dashed', label='SpatialDE', color=colors[0]),
]
lines = lines[::-1]
ll = Legend(ax, lines, [l.get_label() for l in lines], loc='upper left',
            bbox_to_anchor=(1, 1), title='Method')
ax.add_artist(ll)

# plt.show()
plt.savefig('ModuleTiming.svg')

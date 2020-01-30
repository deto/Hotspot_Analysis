import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'

# %% Load Data

data_hs = json.load(open("IDR_hotspot.json"))

data_sde = json.load(open("IDR_SDE.json"))

# %% all pairs in one plot

data_hs = json.load(open("IDR_hotspot.json"))
data_sde = json.load(open("IDR_SDE.json"))

colors = sns.color_palette("deep")
hs_color = colors[0]
sde_color = colors[2]

plt.figure(figsize=(5, 5))

for ds in data_hs:
    x = sorted(np.log10(ds['IDR'])*-1)
    y = len(x) - np.arange(len(x))+1
    handle = plt.plot(x, y, '--', color=hs_color)

for ds in data_sde:
    x = sorted(np.log10(ds['IDR'])*-1)
    y = len(x) - np.arange(len(x))+1
    handle1 = plt.plot(x, y, ':', color=sde_color)

plt.xlim(0, 4)
plt.ylim(0, 1000)
plt.xlabel('IDR ($-log_{10}$)')
plt.ylabel('Genes Above IDR')
plt.legend([handle[0], handle1[0]], ['Hotspot', 'SpatialDE'])
plt.grid(color='#dddddd', dashes=(6, 6))
plt.savefig('IDR_pairs.svg')
# plt.show()

# %%

data_multi = json.load(open("IDR_multi.json"))

plt.figure()

label_map = {
    'hs4': 'Hotspot',
    'hs3': 'Hotspot (without Puck 9)',
    'sd4': 'SpatialDE',
    'sd3': 'SpatialDE (without Puck 9)',
}

for name, ds in data_multi.items():
    x = sorted(np.log10(ds['IDR'])*-1)
    y = len(x) - np.arange(len(x))+1
    plt.plot(x, y, '-', label=label_map[name])

plt.xlim(0, 4)
plt.ylim(0, 1000)
plt.legend()
plt.xlabel('IDR ($-log_{10}$)')
plt.ylabel('Genes Above IDR')
plt.grid(color='#dddddd', dashes=(6, 6))
plt.savefig('IDR_Multi.svg')
#plt.show()

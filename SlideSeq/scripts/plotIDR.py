import json
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'

# %%

data_hs = json.load(open("IDR_hotspot.json"))

plt.figure()

for ds in data_hs:
    x = sorted(np.log10(ds['IDR'])*-1)
    y = len(x) - np.arange(len(x))+1
    plt.plot(x, y, '-')

plt.xlim(0, 8)
plt.ylim(0, 1000)
plt.show()

# %%

data_sde = json.load(open("IDR_SDE.json"))

plt.figure()

for ds in data_sde:
    x = sorted(np.log10(ds['IDR'])*-1)
    y = len(x) - np.arange(len(x))+1
    plt.plot(x, y, '-')

plt.xlim(0, 4)
plt.ylim(0, 1000)
plt.show()

# %% all pairs in one plot

data_hs = json.load(open("IDR_hotspot.json"))
data_sde = json.load(open("IDR_SDE.json"))

colors = plt.get_cmap('tab10').colors

plt.figure()

for ds in data_hs:
    x = sorted(np.log10(ds['IDR'])*-1)
    y = len(x) - np.arange(len(x))+1
    handle = plt.plot(x, y, '--', color=colors[0])

for ds in data_sde:
    x = sorted(np.log10(ds['IDR'])*-1)
    y = len(x) - np.arange(len(x))+1
    handle1 = plt.plot(x, y, ':', color=colors[2])

plt.xlim(0, 4)
plt.ylim(0, 1000)
plt.xlabel('IDR ($-log_{10}$)')
plt.ylabel('Genes Above IDR')
plt.legend([handle[0], handle1[0]], ['Hotspot', 'SpatialDE'])
plt.grid(color='#dddddd', dashes=(6, 6))
plt.savefig('IDR_pairs.svg')

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

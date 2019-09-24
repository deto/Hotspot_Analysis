import loompy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ds = loompy.connect("../../data/SlideSeq/Puck_180819_12/data.loom", mode="r")

scaled = ds.layer['scaled'][:, :]
counts = ds.layer[''][:, :]
symbol = ds.row_attrs['Symbol'][:]
scaled = scaled / 45 * 10000

ds.close()

# %% Calculate stats


mu = scaled.mean(axis=1)
var = scaled.var(axis=1)
disp = (var - mu) / (mu**2)
fano = var / mu
cv2 = var / mu**2

mu_c = counts.mean(axis=1)
var_c = counts.var(axis=1)
disp_c = (var_c - mu_c) / (mu_c**2)
fano_c = var_c / mu_c


# %% Plot

plt.figure()
plt.plot(mu, fano, 'o', ms=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Mean ($\frac{Counts}{10k}$)')
plt.ylabel('Variance / Mean')
plt.show()

# %% Load hotspot results

hs_file = "../Puck_180819_12/hotspot/hotspot.txt"
hs_results = pd.read_table(hs_file, index_col=0)

hs_genes = hs_results.index[hs_results.FDR < .05]
hs_ii = np.array([x in hs_genes for x in symbol])
hs_xx = np.array([x not in hs_results.index for x in symbol])
hs_C_d = {x.Index: x.C for x in hs_results.itertuples()}
hs_C = np.array([hs_C_d[i] if i in hs_C_d else 0 for i in symbol])


# %% Plot w/ HS
colors = plt.get_cmap('tab10').colors

plt.figure(figsize=(4, 4))
plt.plot(mu[hs_xx], fano[hs_xx], 'o', ms=1, color='#DDDDDD')
plt.plot(mu[~hs_xx], fano[~hs_xx], 'o', ms=1, color='#CCCCCC')
ll = plt.plot(mu[hs_ii], fano[hs_ii], 'o', ms=1, color=colors[2])[0]
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Mean ($\frac{Counts}{10k}$)')
plt.ylabel('Variance / Mean')
plt.legend([ll], ['HS-Selected'], markerscale=3, loc='lower right')
plt.subplots_adjust(left=.2, bottom=.2)
# plt.show()
plt.savefig('slideseq_meanvar.svg')

# %% Hover plot

from bio_utils.plots import hover_plot

hover_plot(
    mu,
    var/mu,
    symbol,
    'o', ms=1)
plt.xscale('log')
plt.xlabel(r'Mean ($\frac{Counts}{10k}$)')
plt.ylabel('Variance / Mean')
plt.show()

# %% Plot w/ C

df_all = pd.DataFrame({
    'mu': mu,
    'var': var,
    'fano': var/mu,
    'C': hs_C,
}, index=symbol)

df_all = df_all.sort_values('C')

plt.figure()
plt.scatter(df_all['mu'], df_all['fano'], s=1, c=df_all['C'], vmin=-.05, vmax=.15, cmap="bone_r")
plt.xscale('log')
plt.xlabel(r'Mean ($\frac{Counts}{10k}$)')
plt.ylabel('Variance / Mean')
plt.colorbar()
plt.show()

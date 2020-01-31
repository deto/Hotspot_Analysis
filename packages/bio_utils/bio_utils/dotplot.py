# Going to augment the scatter plot function in seaborn to make a dot plot

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


# Generate some fake data

x_levels = ['A', 'B', 'C']
y_levels = ['D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']


size_data = np.random.uniform(0, 10, size=(len(x_levels), len(y_levels)))
color_data = np.random.uniform(3, 8, size=(len(x_levels), len(y_levels)))


df_rows = []
for i, x in enumerate(x_levels):
    for j, y in enumerate(y_levels):
        df_rows.append([x, y, size_data[i, j], color_data[i, j]])

df = pd.DataFrame(df_rows, columns=['cluster', 'gene', 'zeros', 'exp'])
df = df.sample(df.shape[0])

fig = plt.figure(figsize=(3, 5))

# now, we're in our function
# function arguments

data = df
y = 'gene'
x = 'cluster'
size = 'zeros'
hue = 'exp'
x_order = x_levels
y_order = y_levels
legend_markerscale = 0.5

if fig is None:
    fig = plt.gcf()

if x_order is None:
    x_order = sorted(data[x].unique())
if y_order is None:
    y_order = sorted(data[y].unique())

fig_aspect = fig.get_figheight() / fig.get_figwidth()
data_aspect = len(y_order) / len(x_order)

lpad = 0.1
rpad = 0.05
tpad = 0.05
bpad = 0.1

if data_aspect > fig_aspect:  # height is limiting
    ax_height = 1 - tpad - bpad
    ax_width = ax_height / data_aspect * fig_aspect

else:                         # width is limiting
    ax_width = 1 - lpad - rpad
    ax_height = ax_width * data_aspect / fig_aspect

ax = fig.add_axes([lpad, bpad, ax_width, ax_height])

# Create maps for x, y
x_map = {k: i+.5 for i, k in enumerate(x_order)}
y_map = {k: i+.5 for i, k in enumerate(y_order)}

plot_data = pd.DataFrame({
    x: data[x].map(x_map),
    y: data[y].map(y_map),
    size: data[size],
    hue: data[hue]
})

# Compute maximum circle size
square_side = fig.get_figheight() * ax_height / len(y_order) * 72  # in pts
max_circle_area = (square_side / 2)**2 * 3.14159

with mpl.rc_context(rc={'legend.markerscale': legend_markerscale}):
    sns.scatterplot(
        data=plot_data,
        x=x,
        y=y,
        size=size,
        hue=hue,
        sizes=(0, max_circle_area),
        ax=ax,
        palette='viridis',
    )

plt.xlim(0, len(x_order))
plt.ylim(0, len(y_order))
ax.legend_.set_bbox_to_anchor((1, 1))

ax.set_xticks(np.arange(len(x_order)) + 0.5)
ax.set_xticklabels(x_order)

ax.set_yticks(np.arange(len(y_order)) + 0.5)
ax.set_yticklabels(y_order)

plt.show()

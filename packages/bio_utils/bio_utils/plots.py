import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

__all__ = ["ks_plot", "hover_plot", "dot_plot"]


def ks_plot(*args, **kwargs):
    """
    Arguments should be named Pandas Series
    Keyword argument 'ax' can specify an axes in which to place the plot
    """

    if 'ax' not in kwargs:
        plt.figure()
    else:
        ax = kwargs['ax']
        plt.sca(ax)

    if len(args) == 1 and (isinstance(args, tuple) or isinstance(args, list)):
        args = args[0]

    # Add the minimums/maximums
    min_val = min([min(x) for x in args])
    max_val = max([max(x) for x in args])

    for x1 in args:

        if not isinstance(x1, pd.Series):
            raise ValueError('Inputs should be pandas.Series')

        label1 = x1.name

        x1 = x1.values

        x1 = np.sort(x1)
        y1 = np.arange(0, len(x1)) / len(x1)

        x1 = np.concatenate(([min_val], x1, [max_val]))
        y1 = np.concatenate(([0], y1, [1]))

        plt.plot(x1, y1, label=label1)

    plt.legend(loc='best')


def hover_plot(x, y, labels, *args, ax=None, **kwargs):
    """
    Hover plot with labels in matplotlib

    hover_plot(x, y, labels)

    All other args and kwargs passed into the plt.plot function

    Returns: matplotlib.Figure
    """

    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
    else:
        fig = ax.get_figure()

    line, = plt.plot(x, y, *args, **kwargs)

    annot = ax.annotate("", xy=(0, 0), xytext=(-20, 20),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        x, y = line.get_data()
        annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
        text = "{}".format(" ".join([labels[n] for n in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = line.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

    return fig


def dot_plot(data, x, y, size=None, hue=None,
             square_pad=.9,
             x_order=None, y_order=None, fig=None,
             lpad=0.1, rpad=0.05, tpad=0.05, bpad=0.1,
             **kwargs):

    data = data.copy()

    if fig is None:
        fig = plt.gcf()

    if x_order is None:
        x_order = sorted(data[x].unique())
    if y_order is None:
        y_order = sorted(data[y].unique())
    y_order = y_order[::-1] # Make order from top to bottom

    fig_aspect = fig.get_figheight() / fig.get_figwidth()
    data_aspect = len(y_order) / len(x_order)


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

    data[x] = data[x].map(x_map)
    data[y] = data[y].map(y_map)

    # Compute maximum circle size
    square_side = fig.get_figheight() * ax_height / len(y_order) * 72  # in pts
    max_circle_area = (square_side * square_pad / 2)**2 * 3.14159

    with mpl.rc_context(rc={'legend.markerscale': 1}):
        sns.scatterplot(
            data=data,
            x=x,
            y=y,
            size=size,
            hue=hue,
            sizes=(0, max_circle_area),
            ax=ax,
            **kwargs,
        )

    plt.xlim(0, len(x_order))
    plt.ylim(0, len(y_order))
    ax.legend_.set_bbox_to_anchor((1, 1))

    ax.set_xticks(np.arange(len(x_order)) + 0.5)
    ax.set_xticklabels(x_order)

    ax.set_yticks(np.arange(len(y_order)) + 0.5)
    ax.set_yticklabels(y_order)


# Ok, need two utilities here
# 1) Create colorbar legend
# 2) Create size legend

def create_dot_plot_axes(
    data,
    x,
    y,
    square_pad=0.9,
    x_order=None,
    y_order=None,
    fig=None,
    lpad=0.1,
    rpad=0.05,
    tpad=0.05,
    bpad=0.1,
):

    data = data.copy()

    if fig is None:
        fig = plt.gcf()

    if x_order is None:
        x_order = sorted(data[x].unique())
    if y_order is None:
        y_order = sorted(data[y].unique())

    fig_aspect = fig.get_figheight() / fig.get_figwidth()
    data_aspect = len(y_order) / len(x_order)

    if data_aspect > fig_aspect:  # height is limiting
        ax_height = 1 - tpad - bpad
        ax_width = ax_height / data_aspect * fig_aspect

    else:  # width is limiting
        ax_width = 1 - lpad - rpad
        ax_height = ax_width * data_aspect / fig_aspect

    ax = fig.add_axes([lpad, bpad, ax_width, ax_height])

    # Create maps for x, y
    x_map = {k: i + 0.5 for i, k in enumerate(x_order)}
    y_map = {k: i + 0.5 for i, k in enumerate(y_order)}

    # Compute maximum circle size
    square_side = fig.get_figheight() * ax_height / len(y_order) * 72  # in pts
    max_circle_area = (square_side * square_pad / 2) ** 2 * 3.14159

    plt.xlim(0, len(x_order))
    plt.ylim(0, len(y_order))

    ax.set_xticks(np.arange(len(x_order)) + 0.5)
    ax.set_xticklabels(x_order)

    ax.set_yticks(np.arange(len(y_order)) + 0.5)
    ax.set_yticklabels(y_order)

    return ax, x_map, y_map, max_circle_area



from matplotlib.lines import Line2D

def add_size_legend(sizes, size_to_pts_fn, fig=None, mew=2, **kwargs):
    if fig is None:
        fig = plt.gcf()

    handles = []
    labels = []
    for i in range(len(sizes)):
        labels.append(str(sizes[i]))
        handles.append(
            Line2D([], [], marker='o',
                   markersize=size_to_pts_fn(sizes[i])**.5,
                   markeredgewidth=mew,
                   markeredgecolor='black',
                   markerfacecolor='white',
                   linestyle='None',
                   )
        )

    leg = fig.legend(handles, labels, **kwargs)
    plt.setp(leg.get_title(), multialignment='center')


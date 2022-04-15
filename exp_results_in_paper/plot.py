#! /usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

data_path = # To the unlabeled_time csv
index = ["Dataset", "Query"]
baselines = ["Circinus", "Peregrine", "CBF", "PBE"]
time_limit = 43200

# plot configs
figsize = (32, 7)  # width, height
fontsize = 40
palette = [
    "#293241", "#f28c73", "#537aac", "#81c3d7", "#eb5e28", "#ee6c4d",
    "#e63946", "#e0fbfc"
]
# palette = ["#ff1f5b", "#00cd6c", "#009ade", "#af58ba", "#f28522", "#ffc61e"]
pattern = ['\\', '//', '.', 'x', '-', '\\\\', 'o', '/', 'O', '+']
plt.rc(
    'font', **{
        'size': fontsize,
        'family': 'serif',
        'serif': ['Linux Libertine', 'Times']
    })
plt.rc('text', usetex=True)
plt.rc('axes', prop_cycle=plt.cycler(color=palette))
plt.grid(False)


def plot_function(x, ax):
    # print(df.loc[x])
    plot = df.loc[x].plot.bar(y=baselines,
                              ax=ax,
                              legend=False,
                              edgecolor="white",
                              width=0.8)
    # set filling pattern
    for idx, container in enumerate(ax.containers):
        for bar in container.patches:
            if bar.get_height() == time_limit:
                ax.annotate(r"$\times$",
                            xy=(bar.get_x() + bar.get_width() / 2, time_limit),
                            ha="center",
                            va="bottom",
                            fontsize=fontsize,
                            weight=800)
            bar.set_hatch(pattern[idx % len(pattern)])
    # avoid rotation of x tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    # set frames of subplots
    ax.spines['left'].set_edgecolor(
        'black' if ax.is_first_col() else '#cdcdcd')
    ax.spines['right'].set_edgecolor(
        'black' if ax.is_last_col() else '#cdcdcd')
    ax.yaxis.set_visible(ax.is_first_col())
    ax.set_xlabel(x, weight='bold')
    if ax.is_first_col():
        ax.set_ylabel('Time (s)')
    # ax.set_frame_on(False)
    return plot


def bar_plot(df, grayscale, palette=palette):
    output = "unlabeled_time.pdf"
    if grayscale:
        output = "unlabeled_time_grayscale.pdf"
        palette = [
            "#{0}{0}{0}".format(
                hex(
                    int(
                        int(h[1:3], 16) * 0.299 + int(h[3:5], 16) * 0.587 +
                        int(h[5:7], 16) * 0.114))[2:]) for h in palette
        ]
        print("convert to grayscale", palette)
    plt.rc('axes', prop_cycle=plt.cycler(color=palette))

    n_subplots = len(df.index.levels[0])
    fig, axes = plt.subplots(nrows=1,
                             ncols=n_subplots,
                             sharey=True,
                             figsize=figsize,
                             frameon=False,
                             gridspec_kw={'width_ratios': [2, 2, 2, 1]})
    # for x, ax in zip(df.index.levels[0], axes):
    for x, ax in zip(["YT", "LJ", "OR", "FR"], axes):
        plot_function(x, ax)
    ax = axes[0]
    ax.set_zorder(1)
    ax.legend(prop={'size': fontsize}, loc='upper left', ncol=4, handlelength=2, columnspacing=0.5)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    plt.ylim(0.01, 2e5)
    plt.yscale('log')

    fig.savefig(output, bbox_inches='tight')
    fig.clf()


df = pd.read_csv(data_path).set_index(index)
df = df.clip(upper=time_limit)

bar_plot(df, False)
bar_plot(df, True)

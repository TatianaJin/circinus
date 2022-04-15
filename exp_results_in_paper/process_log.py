#! /usr/bin/env python3
# requires python 3.6+

import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

from argparse import ArgumentParser
from glob import glob
from os import path as osp

plot_configs = {
    "fontsize": 40,
    "palette": ["#f28c73", "#537aac", "#81c3d7", "#eb5e28", "#e0fbfc"],
    "pattern": ['\\', '//', '.', 'x', '|', '\\\\', 'o', '/', 'O', '+']
}

dataset_alias = {
    "human": "HM",
    "youtube": "YT",
    "youtube2007": "YT2",
    "livejournal": "LJ",
    "friendster": "FR",
    "orkut": "OR"
}

dataset_rank = {"human": 0, "youtube": 1, "youtube2007": 2, "livejournal": 3, "orkut": 4, "friendster": 5}


def set_plot_conf():
    plt.rc("font",
           size=plot_configs["fontsize"],
           family="serif",
           serif=['Linux Libertine', 'Times'])
    plt.rc("text", usetex=True)
    plt.rc("axes", prop_cycle=plt.cycler(color=plot_configs["palette"]))
    plt.grid(False)


def get_args():
    parser = ArgumentParser(description="Logs processing")
    subparsers = parser.add_subparsers(dest="command", help="Command")
    # list queries
    list_parser = subparsers.add_parser("list",
                                        help="Get query list from csv logs")
    list_parser.add_argument("logs", nargs="+", help="Paths to the logs.")
    list_parser.add_argument("--em",
                             "--exclude_mismatch",
                             type=int,
                             choices=[0, 1],
                             help="Exclude queries with mismatching results.")
    list_parser.add_argument("-o",
                             "--output",
                             help="Output file name.",
                             default="log")
    # get summary
    summary_parser = subparsers.add_parser(
        "summary", help="Get summary csv from profile logs.")
    summary_parser.add_argument("log_dir", help="Directory of profile logs.")
    summary_parser.add_argument(
        "--mode",
        help="Convert query id to query mode and query index.",
        choices=[1, 0],
        type=int,
        default=True)
    summary_parser.add_argument("-o",
                                "--output",
                                help="Output file name.",
                                default="log")
    summary_parser.add_argument("-i",
                                "--intersection",
                                help="Include intersection count.",
                                choices=[1, 0],
                                type=int,
                                default=True)
    # plot
    plot_parser = subparsers.add_parser("plot", help="Plot summary CSV.")
    plot_parser.add_argument("csv", nargs="+", help="Paths to summary CSVs.")
    plot_parser.add_argument("-m",
                             "--mode",
                             choices=["time", "si", "elapsed", "sum"])
    plot_parser.add_argument("-t",
                             "--titles",
                             nargs="+",
                             help="Titles for summary CSVs.")
    plot_parser.add_argument("-o",
                             "--output",
                             help="Output file name.",
                             default="log")
    plot_parser.add_argument("--legend_pos",
                             help="Position of figure legend",
                             default='lower right')
    plot_parser.add_argument("--ylim", help="Axes ylim", nargs=2, type=float)
    plot_parser.add_argument("--figsize", nargs=2, type=float)
    plot_parser.add_argument("--ratio",
                             help="Compute the ratio based on given data")
    plot_parser.add_argument(
        "-u",
        "--unlabeled",
        action="store_true",
        default=False,
        help="If set, plot for unlabeled queries, else for labeled queries.")
    # plot_uncompressed
    plot_uncompressed_parser = subparsers.add_parser("plot_uncompressed", help="Plot summary CSV for Circinus- and Circinus")
    plot_uncompressed_parser.add_argument("csv", nargs="+", help="Paths to summary CSVs.")
    plot_uncompressed_parser.add_argument("-m", "--mode", choices=["elapsed", "si"], default="elapsed")
    plot_uncompressed_parser.add_argument("--count_threshold", type=int, default=1e8)
    plot_uncompressed_parser.add_argument("-x", choices=["query", "dataset"], default="dataset")
    plot_uncompressed_parser.add_argument("-o", "--output", help="Output file name.", default="log")
    plot_uncompressed_parser.add_argument("--ylim", help="Axes ylim", nargs=2, type=float)
    # plot scaleup
    scaleup_parser = subparsers.add_parser("scaleup", help="Plot scaleup performance")
    scaleup_parser.add_argument("-o", "--output", help="Output file name.", default="log")
    scaleup_parser.add_argument("--legend_pos",
                             help="Position of figure legend",
                             default='lower right')
    # plot motivation
    motivation_parser = subparsers.add_parser("motivation", help="Plot figures for motivation")
    motivation_parser.add_argument("csv", nargs=2, help="Paths to CSVs for unlabeled and labeled queries.")
    motivation_parser.add_argument("-o", "--output", help="Output file name.", default="log")
    motivation_parser.add_argument("--legend_pos", help="Position of figure legend", default='upper right')
    motivation_parser.add_argument("--ylim", help="Axes ylim", nargs=2, type=float, default=(1, 2e13))
    return parser, parser.parse_args()


def list_queries(args):
    queries = None
    for idx, log_path in enumerate(args.logs):
        log = pd.read_csv(log_path).dropna()
        print("{0} complete {1}".format(log_path, len(log)))
        log = log.set_index(
            ["dataset", "query_size", "query_mode", "query_index"])
        log = log[["n_embeddings"]]
        log["n_embeddings"] = pd.to_numeric(log["n_embeddings"])
        if queries is None:
            queries = log
        else:
            queries = queries.join(log, lsuffix="",
                                   rsuffix="_{0}".format(idx)).dropna()
            mismatch = (queries["n_embeddings"] !=
                        queries["n_embeddings_{0}".format(idx)])
            print("{0} mismatches".format(mismatch.sum()))
            if mismatch.sum() > 0:
                print(queries.loc[mismatch])
            if args.em == 1:
                queries = queries.loc[queries["n_embeddings"] == queries[
                    "n_embeddings_{0}".format(idx)]]
    queries = queries[["n_embeddings"]].reset_index()

    queries["id"] = (queries["query_mode"]
                     == "sparse") * 100 + queries["query_index"]
    print(queries.groupby(["dataset", "query_size"]).count())
    queries = queries.set_index(["dataset", "query_size", "id"])
    queries.to_csv(f"{args.output}.csv", columns=["n_embeddings"])


def get_last_line(file_path):
    with open(file_path, 'rb') as f:
        try:  # catch OSError in case of a one line file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        return f.readline().decode().strip()


def get_intersection(file_path):
    intersection = 0.
    with open(file_path, 'r') as f:
        n_ops = int(f.readline().strip().split()[-2])
        for i in range(n_ops - 1):
            line = f.readline()
            if i == 0:
                continue
            intersection += float(line.strip().split(",")[7])
    return intersection


def get_summary(args):
    with open(osp.join(args.log_dir, args.output), "w") as f:
        title = "dataset,query_size,query_mode,query_index,elapsed_execution_time,filter_time,plan_time,enumerate_time,n_embeddings,order,max_task_time"
        if args.intersection:
            f.write(f"{title},n_intersections\n")
        else:
            f.write(f"{title}\n")
        for query in glob(osp.join(args.log_dir, "*.log")):
            log = get_last_line(query)
            if log == "timeout":
                log = ""
            elif args.intersection:
                log = f"{log},{get_intersection(query)}"
            dataset, size, qid = osp.basename(query).split(".")[1:4]
            if args.mode:
                qid = int(qid)
                mode = 'dense' if qid < 100 else 'sparse'
                qid -= (qid >= 100) * 100
                f.write(f"{dataset},{size},{mode},{qid},{log}\n")
            else:
                f.write(f"{dataset},{size},{qid},{log}\n")


def plot_uncompressed(args):
    if args.mode == "elapsed":
        plot_col = 'elapsed_execution_time'
        ylabel = 'Times (s)'
    else:  # args.mode == "si"
        plot_col = 'n_intersections'
        ylabel = 'Intersection Count'
    algo = None
    compared_solution = None
    solutions = ["Circinus-", "Circinus"]
    log_dict = dict({i: None for i in solutions})
    for log_paths in args.csv:
        for log_path in glob(log_paths):
            if "cfl" in log_path:
                algo = "CFL"
            elif "gql" in log_path:
                algo = "GQL"
            else:
                print(f"Cannot get algo of {log_path}, expected CFL or GQL")
            if "none" in log_path and "dynamic" not in log_path:
                compared_solution = "Circinus-"
            else:
                compared_solution = "Circinus"
            log = pd.read_csv(log_path).dropna()
            log = log.loc[log["dataset"] != 'dataset']
            print(f"reading {log_path}, algo={algo}, solution={compared_solution}, complete {len(log)}")
            log["n_embeddings"] = log["n_embeddings"].astype("uint64")
            log = log.loc[log["n_embeddings"] > args.count_threshold]
            log["query_size"] = log["query_size"].astype("uint32")
            log["query_index"] = log["query_index"].astype("uint32")
            log[plot_col] = pd.to_numeric(log[plot_col])
            log["algo"] = algo
            indices = ["dataset", "query_size", "query_mode", "query_index", "algo"]
            log = log.set_index(indices)
            log = log[[plot_col, "n_embeddings"]]
            print(f"{len(log)} samples included")
            log = log.rename(columns={plot_col:compared_solution, "n_embeddings": f"{compared_solution}_count"})
            if log_dict[compared_solution] is None:
                log_dict[compared_solution] = log
            else:
                log_dict[compared_solution] = log_dict[compared_solution].append(log)

    queries = log_dict["Circinus-"].join(log_dict["Circinus"]).dropna()
    queries = queries.loc[queries["Circinus-_count"] == queries["Circinus_count"]]
    queries = queries[solutions]

    queries.columns.name = "solution"
    queries = queries.stack(level=0).reset_index(level=5, name=plot_col)
    if args.x == "dataset":
        queries = queries.reset_index().set_index(["dataset", "algo"])
    else:  # args.x == "query"
        args.x = "query_size"
        queries = queries.reset_index().set_index(["query_size", "algo"])
    print(f"Plot {len(queries)} samples")
    # print(queries)
    set_plot_conf()

    n_subplots = len(queries.index.levels[0])
    print(f"{n_subplots} subplots, {queries.index.levels[1]} x categories")
    figsize = (max(9, len(queries.index.levels[1]) * n_subplots), 8)
    print(f"figsize {figsize}")
    fig, axes = plt.subplots(nrows=1,
                             ncols=n_subplots,
                             sharey=True,
                             figsize=figsize,
                             frameon=False)
    if n_subplots == 1:
        axes = [axes]

    xlabels = queries.index.levels[0]
    if args.x == "dataset":
        xlabels = sorted(xlabels, key=lambda x: dataset_rank[x])
    for x, ax in zip(xlabels, axes):
        q = queries.loc[x].reset_index()
        # print(q)
        box = sns.boxplot(x="algo",
                y=plot_col,
                hue="solution",
                hue_order=solutions,
                data=q,
                ax=ax,
                palette=plot_configs["palette"],
                showfliers=False)
        ax.legend().set_visible(False)
        pattern = plot_configs["pattern"]
        # set filling pattern
        for idx, bar in enumerate(box.patches):  # for legend
            bar.set_hatch(pattern[idx % 2])
        for idx, bar in enumerate(box.artists):
            bar.set_hatch(pattern[idx % 2])
        ax.spines['left'].set_edgecolor(
            'black' if ax.is_first_col() else '#cdcdcd')
        ax.spines['right'].set_edgecolor(
            'black' if ax.is_last_col() else '#cdcdcd')
        ax.yaxis.set_visible(ax.is_first_col())
        ax.set_xlabel(dataset_alias[x] if args.x == "dataset" else x, weight="bold")
        if ax.is_first_col():
            ax.set_ylabel(ylabel)
    # fig.text(-0.075, -0.125, 'Dataset' if args.x == "dataset" else "Query Size", ha='center')
    args.legend_pos = 'lower right'
    if 'left' in args.legend_pos:
        ax = axes[0]
    ax.set_zorder(1)
    lgd = ax.legend(
        prop={'size': plot_configs["fontsize"]},
        handlelength=2,
        labelspacing=0.2,
        columnspacing=0.5,
        loc=args.legend_pos,
        ncol=2)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    if args.ylim is not None:
        plt.ylim(*args.ylim)
    plt.yscale('log')
    print(f"output to {args.output}.pdf")
    fig.savefig(f"{args.output}.pdf",
                bbox_inches='tight',
                bbox_extra_artists=(lgd, ))
    fig.clf()


def plot(args):
    if args.mode == "si":
        args.ylim = (1, 5e9) if args.ylim is None else args.ylim
        plot_col = "n_intersections"
        ylabel = "Intersection Count"
    elif args.mode == "elapsed":
        args.ylim = (5e-4, 1e3) if args.ylim is None else args.ylim
        plot_col = "elapsed_execution_time"
        ylabel = 'Time (s)'
    elif args.mode == "sum":
        args.ylim = (1e-3, 2e2) if "Peregrine" in args.titles else (5e-6, 1e4)
        plot_col = "sum_time"
        ylabel = "Time (s)"
    else:
        args.ylim = (5e-6, 1e4) if args.ylim is None else args.ylim
        plot_col = "enumerate_time"
        ylabel = "Time (s)"
    if args.ratio is not None:
        ylabel = "Ratio over Circinus"
    if args.titles is not None:
        print(args.csv)
        print(args.titles)
        assert len(args.csv) == len(
            args.titles), f"{len(args.csv)}, {len(args.titles)}"
    else:
        args.titles = range(len(args.csv))

    count_threshold = 0 if "Peregrine" in args.titles else 1e8
    queries = None
    log_dict = dict({i: None for i in args.titles})
    for log_title, log_paths in zip(args.titles, args.csv):
        for log_path in glob(log_paths):
            print(f"reading {log_path}")
            log = pd.read_csv(log_path).dropna()
            log = log.loc[log["dataset"] != 'dataset']
            print("{0} complete {1}".format(log_path, len(log)))
            log["query_size"] = log["query_size"].astype("uint32")
            log["query_index"] = log["query_index"].astype("uint32")
            log["n_embeddings"] = log["n_embeddings"].astype("uint64")
            log = log.loc[log["n_embeddings"] > count_threshold]
            if log_title == "Peregrine":
                log["enumerate_time"] = log["enumerate_time"].apply(
                    lambda x: float(x[:-1]))
            log["enumerate_time"] = log["enumerate_time"].astype("float")
            # log = log.loc[log["enumerate_time"] > 1]
            if args.mode == "sum":  # circinus
                if "elapsed_execution_time" in log.columns:
                    log["sum_time"] = log["elapsed_execution_time"] + log[
                        "filter_time"]  # + log["plan_time"]
                elif "build_table_time" in log.columns:  # survey
                    log["sum_time"] = log["nlf_time"] + log[
                        "filter_time"] + log["build_table_time"] + log[
                            "enumerate_time"]  # + log["plan_time"]
                else:  # peregrine
                    log["sum_time"] = log["enumerate_time"]
            log[plot_col] = pd.to_numeric(log[plot_col])
            log = log.set_index(
                ["dataset", "query_size", "query_mode", "query_index"])
            log = log[[plot_col, "n_embeddings"]]
            print(f"{len(log)} samples to include")
            log = log.rename(columns={
                plot_col: f"{log_title}",
                "n_embeddings": f"{log_title}_count"
            })
            if log_dict[log_title] is None:
                log_dict[log_title] = log
            else:
                log_dict[log_title] = log_dict[log_title].append(log)
    previous_title = None
    for title, log in log_dict.items():
        if queries is None:
            queries = log
        else:
            queries = queries.join(log).dropna()
            queries = queries.loc[queries[f"{previous_title}_count"] ==
                                  queries[f"{title}_count"]]
        previous_title = title

    if "Peregrine" in args.titles:
        titles = ["Circinus", "Circinus-CFL", "Peregrine"]
    elif "Circinus-CFL" in args.titles:
        titles = ["Circinus", "Circinus-CFL", "Circinus-GQL", "CFL", "GQL"]
    else:
        titles = sorted(list(set(args.titles)))

    if args.ratio is not None:
        for log_title in titles:
            if log_title != args.ratio:
                queries[log_title] = queries[log_title] / queries[args.ratio]
        titles = [title for title in titles if title != args.ratio]

    queries = queries[set(args.titles)]
    queries.columns.name = "solution"
    queries = queries.stack(level=0).reset_index(level=4, name=plot_col)
    queries = queries.reset_index().set_index(["dataset", "query_size"])
    print(f"Plot {len(queries)} samples")
    # print(queries)

    # pd.set_option('display.max_rows', None)

    set_plot_conf()

    def plot_func(x, ax):
        q = queries.loc[x].reset_index()[["query_size", "solution", plot_col]]
        print(f"{x} records {len(q)}")
        box = sns.boxplot(
            x="query_size",
            y=plot_col,
            hue="solution",
            hue_order=titles,
            data=q,
            ax=ax,
            palette=plot_configs["palette"],
            showfliers=False)  # , medianprops={"color": "white"})
        ax.legend().set_visible(False)
        pattern = plot_configs["pattern"]
        # set filling pattern
        for idx, bar in enumerate(box.patches):  # for legend
            bar.set_hatch(pattern[idx % len(set(args.titles))])
        for idx, bar in enumerate(box.artists):
            bar.set_hatch(pattern[idx % len(set(args.titles))])
        ax.set_xticklabels([
            r"Q$_{{{0}}}$".format(label.get_text())
            for label in ax.get_xticklabels()
        ],
                           rotation=0)
        ax.spines['left'].set_edgecolor(
            'black' if ax.is_first_col() else '#cdcdcd')
        ax.spines['right'].set_edgecolor(
            'black' if ax.is_last_col() else '#cdcdcd')
        ax.yaxis.set_visible(ax.is_first_col())
        ax.set_xlabel(dataset_alias[x], weight="bold")
        if ax.is_first_col():
            ax.set_ylabel(ylabel)

    n_subplots = len(queries.index.levels[0])
    legend_col = 2 if len(titles) == 4 else 3 if len(titles) > 3 else 2
    print(f"{n_subplots} subplots, {queries.index.levels[1]} x categories")
    if args.figsize is not None:
        figsize = args.figsize
    else:
        figsize = (max(8,
                       len(queries.index.levels[1]) * len(titles) * n_subplots /
                       2), 6)
        if "Peregrine" in titles:
            figsize = (14, 6)
            legend_col = 3
    print(f"figsize {figsize}")
    fig, axes = plt.subplots(nrows=1,
                             ncols=n_subplots,
                             sharey=True,
                             figsize=figsize,
                             frameon=False)
    if n_subplots == 1:
        axes = [axes]
    # for x, ax in zip(["youtube2007", "orkut"], axes):
    for x, ax in zip(
            sorted(queries.index.levels[0], key=lambda x: dataset_rank[x]),
            axes):
        plot_func(x, ax)
    if 'left' in args.legend_pos:
        ax = axes[0]
    ax.set_zorder(1)
    lgd = ax.legend(
        prop={'size': plot_configs["fontsize"]},
        frameon=False,
        framealpha=0,
        handlelength=2,
        labelspacing=0.2,
        columnspacing=0.5,
        loc=args.legend_pos,
        ncol=legend_col)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    if args.ylim is not None:
        print(f"ylim {args.ylim}")
        plt.ylim(*args.ylim)
    plt.yscale('log')
    print(f"output to {args.output}.pdf")
    fig.savefig(f"{args.output}.pdf",
                bbox_inches='tight',
                bbox_extra_artists=(lgd, ))
    fig.clf()


def plot_scaleup(args):
    set_plot_conf()
    xticks = [1, 2, 4, 8, 16, 24, 32]
    fig, ax1 = plt.subplots(figsize=(9,9))
    df = pd.read_csv("exp/data/parallelization/memory.csv")
    df["memory"] = df["memory_percentage"] * 128 / 100
    df.plot(x="thread_num", y="memory", ax=ax1, marker='x', label='Peak')
    ax1.set_xlabel("Number of Threads")
    ax1.set_ylabel("Memory Usage (GB)")  # color=
    ax1.set_xticks(xticks)
    ax1.plot(xticks, [22.3158 * 128 / 100]*len(xticks), linestyle='--', label='Graph')
    ax1.set_ylim(0, 40)

    lgd = plt.legend(
            prop={'size': plot_configs["fontsize"]},
            handlelength=2,
            labelspacing=0.2,
            columnspacing=0.5,
            loc=args.legend_pos,
            ncol=3)
    fig.savefig(f"{args.output}_memory.pdf",
                bbox_inches='tight',
                bbox_extra_artists=(lgd, ))
    fig.clf()

    fig, ax1 = plt.subplots(figsize=(9, 9))
    log = None
    for thread_num in xticks:
        path = f"exp/data/parallelization/batch_friendster_16_cfl_cfl_10_dynamic_cc_{thread_num}_parallel_2"
        df = pd.read_csv(path).dropna().set_index(["dataset", "query_size", "query_mode", "query_index"])
        df[thread_num] = df["elapsed_execution_time"]
        df = df[[thread_num]]
        if log is None:
            log = df.loc[df[1] > 30]
        else:
            log = log.join(df).dropna()
            log[thread_num] = log[1] / log[thread_num]

    log[1] = 1
    log.columns.name = "thread_num"
    log = log.stack(level=0).reset_index(level=4, name="speedup")
    print(log)
    sns.boxplot(x="thread_num", y="speedup", data=log, ax=ax1, showfliers=False, color=plot_configs["palette"][0], order=range(1,34))
    ax1.set_xlabel("Number of Threads")
    ax1.set_ylabel("Speedup")
    ax1.set_xticklabels(["1", "2", "", "4", "", "", "", "8", "", "", "", "", "", "", "", "16", "", "", "", "", "", "", "", "24", "", "", "", "", "", "", "", "32", ""])
    #ax1.set_yticks(xticks)
    ax1.set_yticks([2,4,8,16])

    fig.savefig(f"{args.output}_speedup.pdf",
                bbox_inches='tight')
    fig.clf()


def plot_motivation(args):
    set_plot_conf()

    figsize = (8, 8)
    def unlabeled(csv):
        df = pd.read_csv(csv).dropna().set_index(["dataset", "system"])
        df = df.rename(columns={"query1count": "q1", "query3count": "q3"})
        df = df[["q1", "q3"]]
        df.columns.name = "query"
        df = df.stack(level=0).reset_index(level=2, name="si")
        df = df.reset_index()
        df["system"] = df["system"].apply(lambda x: x.capitalize())
        df = df.set_index(["dataset", "query", "system"]).unstack(level=2)
        df.columns = df.columns.droplevel()
        # print(df)

        n_subplots = len(df.index.levels[0])
        fig, axes = plt.subplots(nrows=1,
                                 ncols=n_subplots,
                                 sharey=True,
                                 figsize=figsize,
                                 frameon=False)
        if n_subplots == 1:
            axes = [axes]
        xlabels = sorted(df.index.levels[0], key=lambda x: dataset_rank[x])
        pattern = plot_configs["pattern"]
        for x, ax in zip(xlabels, axes):
            log = df.loc[x].reset_index().set_index(["query"])
            # print(log)
            log.plot.bar(y=["Peregrine", "Circinus", "Optimal"],
                    #width=0.85,
                    ax=ax,
                    legend=False,
                    color=plot_configs["palette"],
                    edgecolor="white")
            ax.set_ylabel("Intersection Count")
            ax.margins(x=0.01, tight=True)
            ax.set_xticklabels([r"$q_{{{0}}}$".format(label.get_text()[-1]) for label in ax.get_xticklabels()], rotation=0)
            ax.spines['left'].set_visible(ax.is_first_col())
            ax.spines['right'].set_visible(ax.is_last_col())
            ax.yaxis.set_visible(ax.is_first_col())
            ax.set_xlabel(dataset_alias[x], weight="bold")
            for idx, container in enumerate(ax.containers):
                for bar in container.patches:
                    bar.set_hatch(pattern[idx % len(pattern)])
        lgd = ax.legend(
            prop={'size': plot_configs["fontsize"]-2},
            handlelength=2,
            labelspacing=0.2,
            columnspacing=0.5,
            loc=args.legend_pos,
            ncol=2)
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        if args.ylim is not None:
            plt.ylim(*args.ylim)
        plt.yscale('log')
        print(f"output to {osp.join(args.output, 'unlabeled_si')}.pdf")
        fig.savefig(f"{osp.join(args.output, 'unlabeled_si')}.pdf",
                    bbox_inches='tight',
                    bbox_extra_artists=(lgd, ))
        fig.clf()

    def labeled(csv):
        df = pd.read_csv(csv).set_index(["dataset"])[["n_intersections", "Mode", "Algorithm"]].dropna()
        df["Mode"] = df["Mode"].apply(lambda x: "Baseline" if x == "UST" else "Circinus" if x == "Circinus" else "Optimal")
        n_subplots = len(df.index.unique())
        fig, axes = plt.subplots(nrows=1,
                                 ncols=n_subplots,
                                 sharey=True,
                                 figsize=figsize,
                                 frameon=False)
        if n_subplots == 1:
            axes = [axes]
        xlabels = sorted(df.index.unique(), key=lambda x: dataset_rank[x])
        pattern = plot_configs["pattern"]
        for x, ax in zip(xlabels, axes):
            log = df.loc[x]
            # print(log)
            box = sns.boxplot(x="Algorithm",
                    y="n_intersections",
                    hue="Mode",
                    hue_order=["Baseline", "Circinus", "Optimal"],
                    ax=ax,
                    showfliers=False,
                    palette=plot_configs["palette"],
                    data=log)
            ax.legend().set_visible(False)
            ax.set_ylabel("Intersection Count")
            ax.margins(x=0.01, tight=True)
            ax.spines['left'].set_visible(ax.is_first_col())
            ax.spines['right'].set_visible(ax.is_last_col())
            ax.yaxis.set_visible(ax.is_first_col())
            ax.set_xlabel(dataset_alias[x], weight="bold")
            for idx, bar in enumerate(box.patches):  # for legend
                bar.set_hatch(pattern[idx % 2])
            for idx, bar in enumerate(box.artists):
                bar.set_hatch(pattern[idx % 2])
        lgd = ax.legend(
            prop={'size': plot_configs["fontsize"]-2},
            handlelength=2,
            labelspacing=0.2,
            columnspacing=0.5,
            loc=args.legend_pos,
            ncol=2)
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        if args.ylim is not None:
            plt.ylim(*args.ylim)
        plt.yscale('log')
        print(f"output to {osp.join(args.output, 'labeled_si')}.pdf")
        fig.savefig(f"{osp.join(args.output, 'labeled_si')}.pdf",
                    bbox_inches='tight',
                    bbox_extra_artists=(lgd, ))
        fig.clf()
    unlabeled(args.csv[0])
    labeled(args.csv[1])


if __name__ == '__main__':
    parser, args = get_args()
    if args.command is None:
        parser.print_help()
        exit(0)
    if args.command == "list":
        list_queries(args)
    elif args.command == "summary":
        get_summary(args)
    elif args.command == "plot":
        plot(args)
    elif args.command == "plot_uncompressed":
        plot_uncompressed(args)
    elif args.command == "scaleup":
        plot_scaleup(args)
    elif args.command == "motivation":
        plot_motivation(args)

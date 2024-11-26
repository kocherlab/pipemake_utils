#!/usr/bin/env python

import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pipemake_utils.misc import confirmFile
from pipemake_utils.logger import startLogger, logArgDict


def argParser():
    parser = argparse.ArgumentParser(description="Filter normalized iHS output")

    parser.add_argument(
        "--input-file",
        help="Input filename",
        type=str,
        action=confirmFile(),
        required=True,
    )
    parser.add_argument(
        "--input-format",
        help="Input file format",
        type=str,
        default="tsv",
        choices=["tsv", "csv"],
    )

    stat_parser = parser.add_mutually_exclusive_group(required=True)
    stat_parser.add_argument("--stat-col", help="Statistic column name", type=str)
    stat_parser.add_argument(
        "--stat-col-int", help="Statistic column location", type=int
    )
    parser.add_argument("--plot-stat-text", help="Statistic plot text", type=str)

    chrom_parser = parser.add_mutually_exclusive_group(required=True)
    chrom_parser.add_argument("--chrom-col", help="CRHOM column name", type=str)
    chrom_parser.add_argument("--chrom-col-int", help="CRHOM column location", type=int)

    pos_parser = parser.add_mutually_exclusive_group(required=True)
    pos_parser.add_argument("--pos-col", help="POS column name", type=str)
    pos_parser.add_argument("--pos-col-int", help="POS column location", type=int)

    parser.add_argument(
        "--chrom-pos-sep",
        help="Chromosome and position separator",
        type=str,
        default="_",
    )

    parser.add_argument("--out-prefix", help="Output prefix", type=str, default="out")

    parser.add_argument(
        "--plot-chrom-text", help="Chromosome plot text", type=str, default="Chromosome"
    )
    parser.add_argument(
        "--plot-abs",
        help="Plot the absolute value of the statistic",
        action="store_true",
    )
    parser.add_argument(
        "--plot-neg-log", help="Plot the -log10 of the statistic", action="store_true"
    )
    parser.add_argument("--plot-dpi", help="Plot DPI", type=int, default=100)

    return vars(parser.parse_args())


def getChromInt(chrom_text, sep):
    # Assign the integer chromosome column location
    chrom_int_pos = None

    for split_pos, split_chrom in enumerate(chrom_text.split(sep)):
        if not split_chrom.isdigit():
            continue

        if chrom_int_pos:
            raise Exception(
                f"Multiple integer chromosome positions found in: {chrom_text}. Using separator: {sep}"
            )

        chrom_int_pos = split_pos

    if not chrom_int_pos:
        raise Exception(
            f"Unable to determine the integer chromosome position from: {chrom_text}. Using separator: {sep}"
        )

    return chrom_int_pos


def main():
    # Parse the arguments
    plot_args = argParser()

    # Start logger and log the arguments
    startLogger(f"{plot_args['out_prefix']}.plot.log")
    logArgDict(plot_args)

    # Assign the input file separator
    if plot_args["input_format"] == "tsv":
        input_separator = "\t"
    elif plot_args["input_format"] == "csv":
        input_separator = ","
    else:
        raise Exception(f"Unknown input format: {plot_args['input_format']}")

    # Read in the data
    if plot_args["chrom_col"] and plot_args["pos_col"]:
        plot_dataframe = pd.read_csv(
            plot_args["input_file"], sep=input_separator, low_memory=False
        )
    else:
        plot_dataframe = pd.read_csv(
            plot_args["input_file"], sep=input_separator, header=None, low_memory=False
        )

    # Assign the plot columns
    plot_chrom_col = (
        plot_args["chrom_col"]
        if plot_args["chrom_col"]
        else plot_args["plot_chrom_text"]
    )
    plot_chrom_int_col = f"{plot_chrom_col}_int"
    plot_pos_col = f"{plot_args['pos_col']}_POS" if plot_args["pos_col"] else "POS"
    plot_stat_col = (
        plot_args["stat_col"] if plot_args["stat_col"] else plot_args["plot_stat_text"]
    )

    # Check if the statistic column needs to be renamed
    if plot_args["stat_col_int"]:
        if not plot_args["plot_stat_text"]:
            raise Exception(
                "--plot-stat-text must be specified if --stat-col-int is specified"
            )
        plot_dataframe = plot_dataframe.rename(
            columns={plot_args["stat_col_int"]: plot_stat_col}
        )

    # Check if the chromosomes and positions are in the same column
    if (
        plot_args["chrom_col_int"] is not None
        and plot_args["chrom_col_int"] == plot_args["pos_col_int"]
    ):
        plot_dataframe[[plot_chrom_col, plot_pos_col]] = plot_dataframe[
            plot_args["chrom_col_int"]
        ].str.rsplit("_", n=1, expand=True)
    elif plot_args["chrom_col"] == plot_args["pos_col"]:
        plot_dataframe[[plot_chrom_col, plot_pos_col]] = plot_dataframe[
            plot_args["chrom_col"]
        ].str.rsplit("_", n=1, expand=True)
    elif plot_args["chrom_col"] and plot_args["pos_col"]:
        plot_dataframe[plot_pos_col] = plot_dataframe[plot_args["pos_col"]].astype(int)

    # Confirm the position column is an integer, for sorting
    plot_dataframe[plot_pos_col] = plot_dataframe[plot_pos_col].astype(int)

    # Reduce the plot dataframe to only the needed columns
    plot_dataframe = plot_dataframe[[plot_chrom_col, plot_pos_col, plot_stat_col]]

    # Plot the absolute value, if specified
    if plot_args["plot_abs"]:
        plot_dataframe[plot_stat_col] = plot_dataframe[plot_stat_col].abs()
        plot_dataframe = plot_dataframe.rename(
            columns={plot_stat_col: f"ABS({plot_stat_col})"}
        )
        plot_stat_col = f"ABS({plot_stat_col})"
    elif plot_args["plot_neg_log"]:
        plot_dataframe[plot_stat_col] = -np.log10(plot_dataframe[plot_stat_col])
        plot_dataframe = plot_dataframe.rename(
            columns={plot_stat_col: f"-log10({plot_stat_col})"}
        )
        plot_stat_col = f"-log10({plot_stat_col})"

    # Check if the chomosome column are integers
    if plot_dataframe[plot_chrom_col].str.isdigit().all():
        plot_dataframe[plot_chrom_int_col] = plot_dataframe[plot_chrom_col].astype(int)
    else:
        # Assign the integer chromosome column location if a separator is specified
        if plot_args["chrom_pos_sep"]:
            chrom_int_pos = getChromInt(
                plot_dataframe[plot_chrom_col][0], plot_args["chrom_pos_sep"]
            )
            plot_dataframe[plot_chrom_int_col] = (
                plot_dataframe[plot_chrom_col]
                .str.split(plot_args["chrom_pos_sep"], expand=True)[chrom_int_pos]
                .astype(int)
            )

        # Otherwise, assign the chromosome column as the integer column
        else:
            plot_dataframe[plot_chrom_int_col] = plot_dataframe[plot_chrom_col]

    # Sort by CHROM and POS
    plot_dataframe = plot_dataframe.sort_values(by=[plot_chrom_int_col, plot_pos_col])

    # Get data length
    plot_data_len = len(plot_dataframe.index)
    plot_dataframe["ORDER"] = range(plot_data_len)

    # Group the dataframe
    grouped_plot_dataframe = plot_dataframe.groupby([plot_chrom_int_col])

    # Create the plot
    plt.rcParams["figure.figsize"] = (24, 6)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colors = ["royalblue", "navy"]
    bg_colors = ["white", "whitesmoke"]
    x_labels = []
    x_labels_pos = []

    for num, (name, group) in enumerate(grouped_plot_dataframe):
        if group.empty:
            continue

        group.plot(
            kind="scatter",
            x="ORDER",
            y=plot_stat_col,
            color=colors[num % len(colors)],
            ax=ax,
            s=1,
            zorder=2,
        )
        plt.axvspan(
            group["ORDER"].iloc[0],
            group["ORDER"].iloc[-1],
            facecolor=bg_colors[num % len(bg_colors)],
            zorder=1,
        )

        x_labels.append(group[plot_chrom_col].iloc[0])
        x_labels_pos.append(
            (
                group["ORDER"].iloc[-1]
                - (group["ORDER"].iloc[-1] - group["ORDER"].iloc[0]) / 2
            )
        )

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, plot_data_len])
    ax.set_ylim(
        [
            plot_dataframe[plot_stat_col].min() * 1.5,
            plot_dataframe[plot_stat_col].max() * 1.5,
        ]
    )
    ax.set_xlabel(plot_chrom_col)

    plt.xticks(rotation=90, fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{plot_args['out_prefix']}.manhattan.png", dpi=plot_args["plot_dpi"])


if __name__ == "__main__":
    main()

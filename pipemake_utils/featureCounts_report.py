#!/usr/bin/env python

import os
import sys

import pandas as pd


def main():
    merged_dataframe = pd.DataFrame()

    for count_file in os.listdir(sys.argv[1]):
        # Filter files
        if count_file.endswith(".summary"):
            continue
        sample = count_file.split(".")[0]

        # Create dataframe and filter
        count_path = os.path.join(sys.argv[1], count_file)
        count_dataframe = pd.read_csv(count_path, sep="\t", comment="#")
        count_dataframe = count_dataframe.set_index("Geneid")

        # Create the series
        count_series = count_dataframe[count_dataframe.columns[-1]]
        count_series = count_series.rename(sample)

        # Merge the series into a dataframe
        if merged_dataframe.empty:
            merged_dataframe = pd.DataFrame(count_series, columns=[sample])
        else:
            merged_dataframe = merged_dataframe.join(count_series, how="outer")

    # Create the file
    merged_dataframe.index.name = None
    merged_dataframe.to_csv(sys.argv[2], sep="\t")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# merge_tables.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import pandas as pd
import glob

def merge_metadata(file):
    df = pd.read_csv(file, dtype=str)  # Ensure all data is read as strings
    df = df.applymap(lambda x: str(x).replace('\r', '').replace('\n', '') if pd.notnull(x) else x)
    return df

# Recursively find all .csv files
csv_files = glob.glob("./**/*.csv", recursive=True)

# Read and merge all CSV files, filling missing columns with NaN
dataframes = [merge_metadata(file) for file in csv_files]
result = pd.concat(dataframes, ignore_index=True, sort=False)

# Write the merged table to a CSV file
result.to_csv("merged.csv", index=False, quoting=1)  # quoting=1 = csv.QUOTE_ALL

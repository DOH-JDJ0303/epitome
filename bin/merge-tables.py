#!/usr/bin/env python3

# merge_tables.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import csv
import glob
import os

def clean_row(row):
    return {k: v.replace('\r', '').replace('\n', '') for k, v in row.items()}

# Recursively find all .csv files
csv_files = glob.glob("./**/*.csv", recursive=True)

# Collect all fieldnames and rows
all_fieldnames = set()
all_rows = []

for file in csv_files:
    with open(file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cleaned = clean_row(row)
            all_fieldnames.update(cleaned.keys())
            all_rows.append(cleaned)

# Write to merged.csv with all combined fieldnames
all_fieldnames = sorted(all_fieldnames)
with open("merged.csv", "w", newline='', encoding='utf-8') as out_f:
    writer = csv.DictWriter(out_f, fieldnames=all_fieldnames, quoting=csv.QUOTE_ALL)
    writer.writeheader()
    for row in all_rows:
        writer.writerow({k: row.get(k, "") for k in all_fieldnames})

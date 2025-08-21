#!/usr/bin/env python
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="Converts the xlsx Tn5 file and checks if the last row has NaNs")
parser.add_argument('-i', '--input', required=True, help="Path to Tn5 file")
parser.add_argument('-o', '--output', required=True, help="Path to save the modified Tn5 file")

args = parser.parse_args()

tn5_barcode_path = args.input

file_extension = os.path.splitext(tn5_barcode_path)[1].lower()

if file_extension == '.csv':
    df = pd.read_csv(tn5_barcode_path)
elif file_extension == '.xlsx' or file_extension == '.xls':
    df = pd.read_excel(tn5_barcode_path)

# Check if the last row has all non-null values (excluding NaN/None)
if df.iloc[-1].isnull().any():
    df = df.iloc[:-1]  # Remove last row if any value is missing

# Save to CSV
df.to_csv(args.output, index=False)



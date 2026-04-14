#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Generate expected demultiplexed samples")
parser.add_argument('-t', '--tn5', required=True, help="Path to Tn5 Barcode Annotation file")
parser.add_argument('-p', '--primer', required=True, help="Path to Primer Annotation file")
parser.add_argument('-o', '--output', required=True, help="Path to save the modified file")

args = parser.parse_args()

tn5_barcode_path = args.tn5
primer_barcode_path = args.primer
output_path = args.output

tn_df = pd.read_csv(tn5_barcode_path)
pb_df = pd.read_csv(primer_barcode_path)
pb_names = pb_df['ID'].drop_duplicates().astype(str).to_numpy()
tn_names = tn_df['Sample Name'].drop_duplicates().astype(str).to_numpy()
cat_names = (tn_names[:, None] + pb_names[None, :]).ravel()
paired = np.array(['_R1.fq.gz', '_R2.fq.gz'])
filename = (cat_names[:, None] + paired[None, :]).ravel()
np.savetxt(output_path, filename, fmt='%s', delimiter='\n')


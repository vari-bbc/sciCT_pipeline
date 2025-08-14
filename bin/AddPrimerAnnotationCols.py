import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Add columns Tn5_s7 and Tn5_s5 to the Primer Annotation index file")
parser.add_argument('-i', '--input', required=True, help="Path to Primer Annotation file")
parser.add_argument('-o', '--output', required=True, help="Path to save the modified file")

args = parser.parse_args()

primer_barcode_path = args.input
pb_df = pd.read_csv(primer_barcode_path)

if "ID" not in pb_df.columns:
    pb_df =df.rename(columns={'Sample': 'ID'})
else: print(f"The column ID is present in the file")

i7_ids = {val: i+1 for i, val in enumerate(pb_df['i7_index_seq'].unique())}
length_i7 = len(i7_ids)

print(f"Number of i7 Ids: {length_i7}" )

pb_df['i7_index_id'] = pb_df['i7_index_seq'].map(i7_ids).apply(lambda x: f"P7_i7_{x}")

i5_ids = {val: i+1 for i, val in enumerate(pb_df['i5_index_seq'].unique())}
length_i5 = len(i5_ids)
print(f"Number of i5 Ids:  {length_i5}" )

pb_df['i5_index_id'] = pb_df['i5_index_seq'].map(i5_ids).apply(lambda x: f"P7_i5_{x}")

pb_df.to_csv(args.output, index= False)
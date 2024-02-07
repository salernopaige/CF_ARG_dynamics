#!/bin/python3
import argparse
import pandas as pd
import os
import re

def normalize_arg_name(arg_name):
    # Remove the _start-stop_+ part from the gene name
    return re.sub(r'_.+$', '', arg_name)

def combine_gcpm(input_files, output_file):
    # Initialize an empty DataFrame for the merged data
    merged_df = pd.DataFrame()

    for file in input_files:
        if os.path.getsize(file) > 0:  # Check if the file is not empty
            try:
                df = pd.read_csv(file, sep='\t')
                if not df.empty and 'ARG_Name' in df.columns:  # Check if DataFrame is not empty and contains 'ARG_Name'
                    # Identify the {sample}_GCPM column
                    gcpm_col = next((col for col in df.columns if col.endswith('_GCPM')), None)
                    if gcpm_col:
                        # Normalize the ARG_Name column
                        df['ARG_Name'] = df['ARG_Name'].apply(normalize_arg_name)

                        # Select only the normalized ARG_Name and {sample}_GCPM columns
                        df = df[['ARG_Name', gcpm_col]]

                        # Perform an outer join on 'ARG_Name'
                        if merged_df.empty:
                            merged_df = df
                        else:
                            merged_df = pd.merge(merged_df, df, on='ARG_Name', how='outer')
                    else:
                        print(f"Warning: No GCPM column found in file {file}. Skipping...")
                else:
                    print(f"Warning: The file {file} is empty or does not contain 'ARG_Name' column. Skipping...")
            except pd.errors.EmptyDataError:
                print(f"Warning: The file {file} is empty and will be skipped.")

    # Group by normalized ARG_Name and sum the GCPM values
    if not merged_df.empty:
        merged_df = merged_df.groupby('ARG_Name', as_index=False).sum()
        merged_df.to_csv(output_file, sep='\t', index=False)
    else:
        print("Warning: No data available to combine. Creating an empty output file.")
        pd.DataFrame().to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine GCPM files into a final table.')
    parser.add_argument('-i', '--input_files', nargs='+', help='List of input GCPM files', required=True)
    parser.add_argument('-o', '--output_file', help='Path to the output combined GCPM file', required=True)
    
    args = parser.parse_args()
    combine_gcpm(args.input_files, args.output_file)

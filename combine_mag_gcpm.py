#!/bin/python3
import argparse
import pandas as pd
import os

def combine_gcpm(input_files, output_file):
    # Initialize an empty DataFrame for the merged data
    merged_df = pd.DataFrame()

    for file in input_files:
        if os.path.getsize(file) > 0:  # Check if the file is not empty
            try:
                df = pd.read_csv(file, sep='\t')
                if not df.empty:  # Check if DataFrame is not empty
                    if 'ARG_Name' in df.columns:
                        # Perform an outer join on 'ARG_Name'
                        merged_df = pd.merge(merged_df, df, on='ARG_Name', how='outer', suffixes=('', '_dup'))
                    else:
                        print(f"Warning: The file {file} does not contain 'ARG_Name' column. Skipping...")
            except pd.errors.EmptyDataError:
                print(f"Warning: The file {file} is empty and will be skipped.")

    # Drop duplicate columns if any
    for col in merged_df.columns:
        if '_dup' in col:
            merged_df.drop(columns=col, inplace=True)

    # Write the combined data to a file
    if not merged_df.empty:
        merged_df.to_csv(output_file, sep='\t', index=False)
    else:
        print("Warning: No data available to combine. Creating an empty output file.")
        pd.DataFrame().to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Combine GCPM files into a final table.')
    parser.add_argument('-i', '--input_files', nargs='+', help='List of input GCPM files', required=True)
    parser.add_argument('-o', '--output_file', help='Path to the output combined GCPM file', required=True)
    
    args = parser.parse_args()
    combine_gcpm(args.input_files, args.output_file)

if __name__ == "__main__":
    main()

#!/bin/python3
import argparse
import pandas as pd

def calculate_mag_gcpm(infile, outfile, sample_name):
    try:
        # Load sample csv
        df = pd.read_csv(infile, sep="\t", names=["ARG_Name", "Sequence_Length", "Mapped_Reads", "Unmapped_Reads"])
    except pd.errors.EmptyDataError:
        print(f"Warning: {infile} is empty. Skipping...")
        return 

    # Ensure numeric data is treated as such
    df['Mapped_Reads'] = pd.to_numeric(df['Mapped_Reads'], errors='coerce')
    df['Sequence_Length'] = pd.to_numeric(df['Sequence_Length'], errors='coerce')
    
    # Filtering out rows where ARG_Name is '*' or Sequence_Length is 0 or data is NA
    df = df.loc[
        (df['ARG_Name'] != '*') & 
        (df['Sequence_Length'] != 0) & 
        df['Mapped_Reads'].notna() & 
        df['Sequence_Length'].notna()
    ]
    
    # If dataframe is empty after filtering, return None
    if df.empty:
        print(f"Warning: No usable data in {infile}. Skipping...")
        return 
    else:
        # Calculating GCPM
        total_ratio = (df['Mapped_Reads'] / df['Sequence_Length']).sum()
        df['GCPM'] = ((df['Mapped_Reads'] / df['Sequence_Length']) * 1e6) / total_ratio

        # Assertion to check if the sum of GCPM values is close to 1e6
        gcpm_sum = df['GCPM'].sum()
        assert round(gcpm_sum, 0) == 1e6, f"Sum of GCPM values for {sample_name} is {gcpm_sum}, not equal to 1e6"
        
        df.rename(columns={'GCPM': f'{sample_name}_GCPM'}, inplace=True)
        df.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate GCPM values for each gene in each sample, then iteratively combine GCPM values from each sample to final table')
    parser.add_argument('-i', '--input', help='Path to input files', required=True)
    parser.add_argument('-o', '--output', help='Path to output file', required=True)
    parser.add_argument('-s', '--sample_name', help = 'Sample Name', required = True)
    args = parser.parse_args()
    
    calculate_mag_gcpm(args.input, args.output, args.sample_name)


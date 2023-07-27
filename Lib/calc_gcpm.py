import argparse
import pandas as pd
import numpy as np

def calculate_gcpm(counts_table, lengths_table):
    lengths_table['ARO_Term'] = lengths_table['ARO_Term'].str.strip()
    counts_table['ARO Accession'] = counts_table['ARO Accession'].astype(str)
    lengths_table['ARO_Term'] = lengths_table['ARO_Term'].astype(str)
    
    # Copy the counts_table to refill in later
    gcpm_table=counts_table.copy()

    # Divide all read counts by respective gene lengths
    merged = pd.merge(gcpm_table, lengths_table, left_on='ARO Accession', right_on='ARO_Term')
    total_counts = merged[['ARO_Term']].copy()
    
    # ID columns with count information 
    count_columns = [col for col in counts_table.columns if col.startswith("SRR") and col.endswith("_Mapped_Reads")]

    # Calculate total counts for each sample
    for sample in count_columns:
        new_col_name = sample.split('_')[0]
        total_counts.loc[:, new_col_name] = merged[sample] / merged['Gene_Length']

    count_info = {}
    aro_accessions_count_info = set()  # To store unique ARO accessions from count_info

    for gene_id, gene_row in counts_table.iterrows():
        aro_accession = gene_row["ARO Accession"]
        for sample in count_columns:
            count_value = gene_row[sample]
            if sample not in count_info:
                count_info[sample] = {'ARO Accession': [], 'Count Value': []}
            count_info[sample]['ARO Accession'].append(aro_accession)
            count_info[sample]['Count Value'].append(count_value)
            aro_accessions_count_info.add(aro_accession)  # Add the ARO accession to the set

    gcpm_values = {}

    for sample, sample_data in count_info.items():
        aro_accessions = sample_data['ARO Accession']
        counts = sample_data['Count Value']
        sample_total_counts = total_counts[sample.split('_')[0]].sum()
        new_col = sample.split('_')[0] + '_GCPM'
        gcpm_table.rename({sample:new_col}, inplace=True)

        for aro_accession, count in zip(aro_accessions, counts):
            gene_ARO = lengths_table[lengths_table['ARO_Term'] == aro_accession]
            if not gene_ARO.empty:    
                gene_length = gene_ARO['Gene_Length'].values[0]
            else:
                print(f"No matching ARO accession in lengths_table: {aro_accession}")

            if count > 0:
                gcpm = ((count / gene_length) * 10**6) / (sample_total_counts)
                # Replace the column in gcpm_table with the GCPM values
                row = gcpm_table[gcpm_table['ARO Accession'] == aro_accession].index[0]
                gcpm_table.loc[row, new_col] = gcpm
                # print(f"The GCPM for {aro_accession} is {gcpm}")
            else:
                gcpm = 0
                # Replace the column in gcpm_table with the GCPM values
                row = gcpm_table[gcpm_table['ARO Accession'] == aro_accession].index[0]
                gcpm_table.loc[row, new_col] = gcpm
                # print(f"Gene count for {aro_accession} is 0, no GCPM can be calculated")

        gcpm_table.drop(columns=sample, inplace=True)

    return(gcpm_table)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate GCPM values for each sample and each gene in a Mapped Reads file')
    parser.add_argument('-c', '--counts_table', help='Path for the Mapped Reads file', required=True)
    parser.add_argument('-g', '--gene_lengths', help='Path for the Gene Lengths file', required=True)
    parser.add_argument('-o', '--outfile', help='Path for output file', required=True)
    args = parser.parse_args()

    # Read the input tables
    counts_table = pd.read_csv(args.counts_table).fillna(0)
    count_columns = [col for col in counts_table.columns if col.startswith("SRR") and col.endswith("_Mapped_Reads")]
    counts_table[count_columns] = counts_table[count_columns].apply(pd.to_numeric, errors='coerce')  # Convert numeric columns to numeric, non-numeric columns will remain as objects
    gene_lengths_table = pd.read_csv(args.gene_lengths)

    # Remove "ARO:" prefix from "ARO_Term" column in gene_lengths_table
    gene_lengths_table["ARO_Term"] = gene_lengths_table["ARO_Term"].str.replace("ARO:", "")
    gene_lengths_table["Gene_Length"] = pd.to_numeric(gene_lengths_table["Gene_Length"], errors='coerce')  # Convert "Gene_Length" column to numeric

    # Calculate GCPM values
    gcpm_table = calculate_gcpm(counts_table, gene_lengths_table)

    # Save the output table
    gcpm_table.to_csv(args.outfile, index=False)

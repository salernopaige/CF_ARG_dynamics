import argparse
import pandas as pd
import numpy as np

def calculate_gcpm(counts_table, lengths_table):	
	# Copy the counts_table to refill in later
	gcpm_table=counts_table.copy()

	# Merge the dataframes
	merged = pd.merge(gcpm_table, lengths_table, left_on='ARO Accession', right_on='ARO_Term')   
	merged.drop(columns=['ARO_Term', 'Gene_Name'], inplace=True)

	# ID columns with count information 
	count_columns = [col for col in counts_table.columns if col.startswith("SRR") and col.endswith("_Mapped_Reads")]

	# Calculate total counts for each sample
	for sample in count_columns:
		# Split sample name for new column names 
		new_col_name = sample.split('_')[0]
		ratio_name = new_col_name + '_ratio'		
		gcpm_name = new_col_name + '_gcpm'

		# Calculate the ratio of counts to gene length for each gene 
		merged[ratio_name] = merged[sample] / merged['Gene_Length']

		# Calculate the sum of ratios for all genes
		sum_ratios = merged[ratio_name].sum()

		# Calculate the 'gcpm' value for each gene
		merged[gcpm_name] = (merged[ratio_name] * 1e6) / sum_ratios
 
        # Store the sum of GCPM values for each sample
		gcpm_sum = merged[gcpm_name].sum()
		gcpm_sum = round(gcpm_sum, 1)

	    # Check if the sum of GCPM values for this sample is equal to 1e6
		assert gcpm_sum == 1e6, f"Sum of GCPM values for {new_col_name} is {gcpm_sum}, not equal to 1e6"
	
		# Drop unwanted columns to get final table
		merged.drop(columns=sample, inplace=True)
		merged.drop(columns=ratio_name, inplace=True)
		print(merged.columns)
	
	return(merged)

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
	counts_table['ARO Accession'] = counts_table['ARO Accession'].astype(str)

	# Remove "ARO:" prefix from "ARO_Term" column in gene_lengths_table
	gene_lengths_table = pd.read_csv(args.gene_lengths)
	gene_lengths_table["ARO_Term"] = gene_lengths_table["ARO_Term"].str.replace("ARO:", "")
	gene_lengths_table["Gene_Length"] = pd.to_numeric(gene_lengths_table["Gene_Length"], errors='coerce')  # Convert "Gene_Length" column to numeric
	gene_lengths_table['ARO_Term'] = gene_lengths_table['ARO_Term'].str.strip()
	gene_lengths_table['ARO_Term'] = gene_lengths_table['ARO_Term'].astype(str)

	# Calculate GCPM values
	gcpm_table = calculate_gcpm(counts_table, gene_lengths_table)

	# Save the output table
	gcpm_table.to_csv(args.outfile, index=False)
